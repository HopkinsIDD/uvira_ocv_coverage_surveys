## code to run stan model 

# Stan model assumes the coverage rate decreases over time exponentially with a rate that is proportional to the current coverage rate

# Load libraries
library(cmdstanr)
library(bayesplot)
library(here)
library(tidyverse)
library(here)

#https://mc-stan.org/cmdstanr/
# Need to run these two steps:
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#cmdstanr::install_cmdstan()

# Andrew version
file <- file.path(here("notebooks/uvira_pop_movement/population_movement_v3.stan"))

# Aybüke version
#file <- file.path(here("notebooks/uvira_pop_movement/population_movement_v3_old.stan"))
mod <- cmdstan_model(file) 

# Define the number of age groups
num_age_groups <- 6 # decision to remove the <1 category on Dec 23
num_groups <- num_age_groups
num_surveys <- 3  # Number of survey rounds

# Define the age and sex categories
age_categories <- 1:num_age_groups

# Define the time values for each survey round
time_round1 <- 11
time_round2 <- 19 # updated in minor revision October 2024
time_round3 <- 32 # updated in minor revision October 2024

# Construct the t variable
t <- matrix(rep(c(time_round1, time_round2, time_round3), each = num_groups), nrow = num_surveys, ncol = num_groups, byrow = TRUE)

# Define the coverage rates for each age group and survey round
# Two or more doses of OCV
coverage <- matrix(c(
  # Round 1: Coverage for each age group
 0.26, 0.27, 0.21, 0.25, 0.21, 0.21,
  # Round 2: Coverage for each age group
 0.23, 0.24, 0.20, 0.20, 0.22, 0.18,
    # Round 3: Coverage for each age group
 0.11, 0.10, 0.08, 0.10, 0.16, 0.19), nrow=num_surveys, ncol=num_groups, byrow=TRUE)


# Define the population sizes for each age group based 
# on census data in that survey round
population_sizes <- matrix(c(
  # Round 1: 6 age groups. Census from 2021
  (14652119), (16294305), (24554072), (24015985), (12739523), (8777814),
  # Round 2: 6 age groups. Census from 2022
  (15017237), (16609090), (25589879), (24764164), (13237658), (9092663),
  # Round 3: 6 age groups. Census from 2023
  (15378752), (16965832), (26627701), (25519326), (13755051), (9427524)), nrow=num_surveys, ncol=num_groups, byrow=TRUE)


# updated in minor revision October 2024: updated to only include individuals not missing vax status
# Define the number of sampled individuals for each age group and survey round
sampled_individuals <- matrix(c(
  # Round 1
  346, 403, 641, 437, 210, 179,
  # Round 2
  554, 578, 967, 636, 325, 246,
  # Round 3
  432, 500, 745, 430, 241, 157), nrow = 3, ncol = num_groups, byrow = TRUE)

# sampled_individuals <- matrix(c(
#   # Round 1
#   354, 405, 642, 438, 210, 180,
#   # Round 2
#   569, 585, 972, 641, 326, 250,
#   # Round 3
#   439, 514, 759, 446, 244, 160), nrow = 3, ncol = num_groups, byrow = TRUE)
# 


# Compute the number vaccinated for each age/sex group and survey round
num_vaccinated <- coverage*sampled_individuals

# Create data_list with age, sex, and population size information
data_list <- list(
  S = num_surveys,  # Number of survey rounds
  G = num_groups,  # Number of age/sex groups
  v = num_vaccinated ,  # Number vaccinated in each age group
  n = sampled_individuals,  # Number of sampled individuals for each age/sex group and survey round
  t = t,  # Time since vaccination for each agex group and survey round
  N = population_sizes  # Population size for each age group
)


# Run the Stan model
# Output should be the predicted number vaccinated in each group for each time point 
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500  # Print update every 500 iterations
)

draws <- fit$draws(format = "df")

mcmc_hist(fit$draws("lambda"))

predicted_coverage <- fit$draws("coverage", format = "df")

#make into long format and parse the time from the coverage column names
predicted_coverage <- predicted_coverage %>%
  pivot_longer(cols = starts_with("coverage"), names_to = "group", values_to = "coverage")

# pull out time and group variables
predicted_coverage <- predicted_coverage %>%
  mutate(
    time = as.numeric(str_extract(group, "(?<=\\[)[^,]+")),
    group = as.integer(str_extract(group, "(?<=,)[^\\]]+"))
  )

quantiles_per_time <- predicted_coverage |>
  group_by(time, group) |>
  summarise(q2_5 = quantile(coverage, probs = 0.025, na.rm=TRUE),
            q5 = quantile(coverage, probs = 0.05, na.rm=TRUE),
            q25 = quantile(coverage, probs = 0.25, na.rm=TRUE),
            q50 = quantile(coverage, probs = 0.5, na.rm=TRUE),
            q75 = quantile(coverage, probs = 0.75, na.rm=TRUE),
            q95 = quantile(coverage, probs = 0.95, na.rm=TRUE),
            q975 = quantile(coverage, probs = 0.975, na.rm=TRUE))

quantiles_per_time$group <- factor(quantiles_per_time$group)

# Create an overall coverage 
# Define the population sizes for each age group. Overall, across all 3 survey rounds
population_sizes <- c(1332, 1481, 2353, 1503, 776, 582) # corrected

# # Calculate the weighted average coverage decline for each time point
coverage_overall <- quantiles_per_time %>%
  group_by(time) %>%
  summarise(weighted_average = sum(q50 * population_sizes) / sum(population_sizes))


# Try merging back in the weighted average for each time point 
quantiles_per_time <- quantiles_per_time %>% inner_join(coverage_overall, by = c("time" = "time"))

# Calculate standard error
quantiles_per_time <- quantiles_per_time %>%
  mutate(
    standard_error = sqrt(
      sum((q50 - weighted_average)^2 * population_sizes) /
        (sum(population_sizes) * (length(population_sizes) - 1))))


# Calculate t-score for the desired confidence level (e.g., 95%)
confidence_level <- 0.95
t_score = qt((1 + confidence_level) / 2, df = length(population_sizes) - 1)

# Calculate margin of error
quantiles_per_time <- quantiles_per_time %>%
  mutate(
    margin_of_error = t_score * standard_error)

# Calculate lower and upper bounds of the confidence interval
quantiles_per_time <- quantiles_per_time %>%
  mutate(
    lower_ci = weighted_average - margin_of_error,
    upper_ci = weighted_average + margin_of_error)

# Adjust the indexing
quantiles_per_time <- quantiles_per_time %>% mutate(time=time-1)

coverage_overall_final <- quantiles_per_time %>% distinct(time, weighted_average, .keep_all = TRUE) %>% select(time, weighted_average, lower_ci, upper_ci)

# # Plot overall decline with weighted average
# plot_overall <- ggplot(data = coverage_overall_final, aes(x = time, y = weighted_average)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "gray", linewidth=0.1) +
#   geom_line(aes(x = time, y = weighted_average)) +
#   labs(x = "Months post-vaccination", y = "At least 2-dose coverage", title = "") +
#   geom_vline(xintercept = data_list$t[1], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[2], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[3], linetype = "dashed") + theme_minimal() + ylim(0, 1.0)
# 
# 
# plot_overall
# 
# # Export
# ggsave("decline_stan_overall_2dose.pdf", plot = plot_overall, width = 10, height = 7)
# 

# # Plot age specific decline
# plot_age <- ggplot(data = quantiles_per_time, aes(x = time, y = q50, color = group)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = q2_5, ymax = q975), alpha = 0.3, fill = "gray", linewidth=0.1) +
#   geom_line(aes(x = time, y = q50)) +
#   labs(x = "Months post-vaccination", y = "At least 2-dose coverage", title = "") +
#   geom_vline(xintercept = data_list$t[1], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[2], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[3], linetype = "dashed") +
#   scale_color_discrete(name = "Age at time of campaign", labels = c("1-4", "5-9", "10-19", "20-34", "35-49", "50+")) + 
#   theme_minimal() + ylim(0, 1.0) + theme(legend.position = "bottom")
# 
# 
# plot_age
# ggsave("decline_stan_age_2dose.pdf", plot = plot_age, width = 10, height = 7)

# Create overall coverage points
# Remove those <1 
coverage_master <- master  %>% filter(is.na(n_dose2plus)==FALSE) %>% filter(age_cat_elig3!="-1.762-0") %>% 
  filter(is.na(age_cat_elig3)==FALSE) %>% group_by(time_since) %>%
  summarize(mean_coverage= mean(n_dose2plus))

# Create age-specific coverage points
coverage_master_age <- master %>% filter(is.na(n_dose2plus)==FALSE) %>% filter(age_cat_elig3!="-1.762-0") %>% 
  filter(is.na(age_cat_elig3)==FALSE) %>% group_by(age_cat_elig3, time_since) %>%
  summarize(mean_coverage= mean(n_dose2plus)) %>% mutate(group=case_when(
    age_cat_elig3=="1-4" ~ 1,
    age_cat_elig3=="5-9" ~ 2,
    age_cat_elig3=="10-19" ~ 3,
    age_cat_elig3=="20-34" ~ 4,
    age_cat_elig3=="35-49" ~ 5,
    age_cat_elig3=="50+" ~ 6,
    TRUE ~ NA_real_))

coverage_master_age$group <- factor(coverage_master_age$group)

# Overall with points

plot_points <- ggplot(data = coverage_overall_final, aes(x = time, y = weighted_average)) +
  geom_line() +
  geom_point(data = coverage_master, aes(x = time_since, y = mean_coverage), size=3, alpha=0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "gray", linewidth=0.1) +
  geom_line(aes(x = time, y = weighted_average)) +
  labs(x = "Months post-vaccination", y = "At least 2-dose coverage", title = "") +
  geom_vline(xintercept = data_list$t[1], linetype = "dashed") +
  geom_vline(xintercept = data_list$t[2], linetype = "dashed") +
  geom_vline(xintercept = data_list$t[3], linetype = "dashed") + theme_minimal() + ylim(0, 1.0)

plot_points

# Export
ggsave("decline_stan_points_2dose.pdf", plot = plot_points, width = 16, height = 7)

# # Stratified by age with points
# plot_points_age <- ggplot(data = quantiles_per_time, aes(x = time, y = q50, color = group)) +
#   geom_ribbon(aes(ymin = q2_5, ymax = q975), alpha = 0.3, fill = "gray", linewidth=0.1) +
#   geom_line(aes(x = time, y = q50)) +
#   labs(x = "Months post-vaccination", y = "At least 2-dose coverage", title = "") +
#   geom_vline(xintercept = data_list$t[1], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[2], linetype = "dashed") +
#   geom_vline(xintercept = data_list$t[3], linetype = "dashed") +
#   scale_color_discrete(name = "Age at time of campaign", labels = c("1-4", "5-9", "10-19", "20-34", "35-49", "50+")) +
#   geom_point(data = coverage_master_age, aes(x = time_since, y = mean_coverage, color = group), size=3, alpha=0.7) +
#   theme_minimal() + theme(legend.position = "bottom") + ylim(0, 1.0)
# 
# plot_points_age
# 
# # Export
# ggsave("decline_stan_age_points.pdf", plot = plot_points_age, width = 10, height = 7)

# Try scaling points based on population size
coverage_master_age <- coverage_master_age %>% mutate(pop_size = case_when(
  time_since==11 & group==1 ~ 346, 
  time_since==11 & group==2 ~ 403, 
  time_since==11 & group==3 ~ 641, 
  time_since==11 & group==4 ~ 437, 
  time_since==11 & group==5 ~ 210, 
  time_since==11 & group==6 ~ 179, 
  time_since==19 & group==1 ~ 554,
  time_since==19 & group==2 ~ 578,
  time_since==19 & group==3 ~ 967,
  time_since==19 & group==4 ~ 636,
  time_since==19 & group==5 ~ 325,
  time_since==19 & group==6 ~ 246,
  time_since==32 & group==1 ~ 432,
  time_since==32 & group==2 ~ 500,
  time_since==32 & group==3 ~ 745,
  time_since==32 & group==4 ~ 430,
  time_since==32 & group==5 ~ 241,
  time_since==32 & group==6 ~ 157,
  TRUE ~ NA ))


plot_points_age_scaled <- ggplot(data = quantiles_per_time, aes(x = time, y = q50, color = group)) +
  geom_ribbon(aes(ymin = q2_5, ymax = q975), alpha = 0.3, fill = "gray", linewidth=0.1) +
  geom_line(aes(x = time, y = q50)) +
  labs(x = "Months post-vaccination", y = "At least 2-dose coverage", title = "") +
  geom_vline(xintercept = data_list$t[1], linetype = "dashed") +
  geom_vline(xintercept = data_list$t[2], linetype = "dashed") +
  geom_vline(xintercept = data_list$t[3], linetype = "dashed") +
  scale_color_discrete(name = "Age at time of campaign", labels = c("1-4", "5-9", "10-19", "20-34", "35-49", "50+")) +
  scale_size_continuous(name = "Sampled population size") +  
  geom_point(data = coverage_master_age, aes(x = time_since, y = mean_coverage, size = pop_size)) +
  geom_point(data = coverage_master_age, aes(x = time_since, y = mean_coverage, color = group, size = pop_size), alpha=0.7) +
  theme_minimal() + theme(legend.position = "bottom",
                          legend.key.size = unit(0.3, "lines"),  # Adjust key size
                          legend.text = element_text(size = 5),  # Adjust text size
                          legend.title = element_text(size = 6)  # Adjust title size 
  ) + ylim(0, 1.0) + guides(color = guide_legend(order = 2), size = guide_legend(order = 1))


plot_points_age_scaled

# Put the plots together

library(patchwork)
#full <- plot_points + plot_points_age
full <- plot_points + plot_points_age_scaled
ggsave("decline_stan_overall_age_2dose.pdf", plot = full, width = 16, height = 7)


# Determine average decay rate by analyzing posterior distribution of the estimated lambda parameter
# Calculate summary statistics
median_lambda_1 <- median(draws$'lambda[1]')
mean_lambda_1 <- mean(draws$'lambda[1]')
credible_interval_lambda_1 <- quantile(draws$'lambda[1]', c(0.025, 0.975))

median_lambda_2 <- median(draws$'lambda[2]')
mean_lambda_2 <- mean(draws$'lambda[2]')
credible_interval_lambda_2 <- quantile(draws$'lambda[2]', c(0.025, 0.975))

median_lambda_3 <- median(draws$'lambda[3]')
mean_lambda_3 <- mean(draws$'lambda[3]')
credible_interval_lambda_3 <- quantile(draws$'lambda[3]', c(0.025, 0.975))

median_lambda_4 <- median(draws$'lambda[4]')
mean_lambda_4 <- mean(draws$'lambda[4]')
credible_interval_lambda_4 <- quantile(draws$'lambda[4]', c(0.025, 0.975))

median_lambda_5 <- median(draws$'lambda[5]')
mean_lambda_5 <- mean(draws$'lambda[5]')
credible_interval_lambda_5 <- quantile(draws$'lambda[5]', c(0.025, 0.975))

median_lambda_6 <- median(draws$'lambda[6]')
mean_lambda_6 <- mean(draws$'lambda[6]')
credible_interval_lambda_6 <- quantile(draws$'lambda[6]', c(0.025, 0.975))

# Print summary statistics
cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_1[1], "-", credible_interval_lambda_1[2], "\n")

cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_2[1], "-", credible_interval_lambda_2[2], "\n")

cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_3[1], "-", credible_interval_lambda_3[2], "\n")

cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_4[1], "-", credible_interval_lambda_4[2], "\n")

cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_5[1], "-", credible_interval_lambda_5[2], "\n")

cat("95% Credible Interval for Decay Rate (lambda):", 
    credible_interval_lambda_6[1], "-", credible_interval_lambda_6[2], "\n")


# Estimate yearly decline instead of monthly decline
mean_lambda_1 = mean_lambda_1*12
mean_lambda_2 = mean_lambda_2*12
mean_lambda_3 = mean_lambda_3*12
mean_lambda_4 = mean_lambda_4*12
mean_lambda_5 = mean_lambda_5*12
mean_lambda_6 = mean_lambda_6*12


# Given confidence interval endpoints in months
lower_ci_per_month1 <- credible_interval_lambda_1[1]
upper_ci_per_month1 <- credible_interval_lambda_1[2]

lower_ci_per_month2 <- credible_interval_lambda_2[1]
upper_ci_per_month2 <- credible_interval_lambda_2[2]

lower_ci_per_month3 <- credible_interval_lambda_3[1]
upper_ci_per_month3 <- credible_interval_lambda_3[2]

lower_ci_per_month4 <- credible_interval_lambda_4[1]
upper_ci_per_month4 <- credible_interval_lambda_4[2]

lower_ci_per_month5 <- credible_interval_lambda_5[1]
upper_ci_per_month5 <- credible_interval_lambda_5[2]

lower_ci_per_month6 <- credible_interval_lambda_6[1]
upper_ci_per_month6 <- credible_interval_lambda_6[2]


# Convert to annual decay rate confidence interval
lower_ci_per_year1 <- lower_ci_per_month1 * 12
upper_ci_per_year1 <- upper_ci_per_month1 * 12

lower_ci_per_year2 <- lower_ci_per_month2 * 12
upper_ci_per_year2 <- upper_ci_per_month2 * 12

lower_ci_per_year3 <- lower_ci_per_month3 * 12
upper_ci_per_year3 <- upper_ci_per_month3 * 12

lower_ci_per_year4 <- lower_ci_per_month4 * 12
upper_ci_per_year4 <- upper_ci_per_month4 * 12

lower_ci_per_year5 <- lower_ci_per_month5 * 12
upper_ci_per_year5 <- upper_ci_per_month5 * 12

lower_ci_per_year6 <- lower_ci_per_month6 * 12
upper_ci_per_year6 <- upper_ci_per_month6 * 12

# Print the result
cat("Mean Decay Rate (lambda):", mean_lambda_1, "\n")
cat("Mean Decay Rate (lambda):", mean_lambda_2, "\n")
cat("Mean Decay Rate (lambda):", mean_lambda_3, "\n")
cat("Mean Decay Rate (lambda):", mean_lambda_4, "\n")
cat("Mean Decay Rate (lambda):", mean_lambda_5, "\n")
cat("Mean Decay Rate (lambda):", mean_lambda_6, "\n")


cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year1, "-", upper_ci_per_year1, "\n")
cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year2, "-", upper_ci_per_year2, "\n")
cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year3, "-", upper_ci_per_year3, "\n")
cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year4, "-", upper_ci_per_year4, "\n")
cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year5, "-", upper_ci_per_year5, "\n")
cat("Annual Decay Rate Confidence Interval:", lower_ci_per_year6, "-", upper_ci_per_year6, "\n")

