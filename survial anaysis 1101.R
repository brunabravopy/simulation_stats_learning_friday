# AFT Weibull Simulation with Multiple Sample Sizes and Censoring Schemes
# Model: T = T_0 * exp(η_AFT) where T_0 ~ Weibull(shape=α, scale=β)
# η_AFT = b_age * age_binary + b_smoker * smoker + b_sex * sex
# Covariates: age (0=<70, 1=≥70), smoker (0=non-smoker, 1=smoker), sex (0=female, 1=male)

# For each covariate combination, generate n_points patients from that distribution
# Then apply 3 censoring times to each dataset

###Generate the data
set.seed(3456)

# Sample sizes to generate
n_points <- c(100, 300, 500)

# Censoring times (study end dates)
time1 <- 5  # years
time2 <- 10  # years  
time3 <- 20  # years
censoring_times <- c(time1, time2, time3)

# AFT Weibull parameters
alpha <- 2.1                 # Weibull shape parameter
L_med <- 15                  # Desired median lifetime for baseline group (years)
beta <- L_med / (log(2))^(1 / alpha)  # Weibull scale parameter

# AFT coefficients (for binary covariates)
b_age <- 0.2                 # Age effect (1=≥70 years vs 0=<70 years)
b_smoker <- -0.3              # Smoking effect (1=smoker vs 0=non-smoker)
b_sex <- 0.1                 # Sex effect (1=male vs 0=female)

# All 8 covariate combinations, all binary covariates
covariate_combos <- expand.grid(age = 0:1, smoker = 0:1, sex = 0:1)

#make names in the covariate_combos dataframe --> used for dataset naming later
covariate_combos$combo_label <- paste0(
  "age", covariate_combos$age, 
  "_smoker", covariate_combos$smoker,
  "_sex", covariate_combos$sex
)

generate_AFT_times <- function(n, age_val, smoker_val, sex_val, 
                               alpha, beta, b_age, b_smoker, b_sex) {
  # Calculate η_AFT for this covariate combination (same for all n patients)
  eta_AFT <- b_age * age_val + b_smoker * smoker_val + b_sex * sex_val
  
  T_0 <- rweibull(n, shape = alpha, scale = beta)
  
  points <- T_0 * exp(eta_AFT)
  
  return(list(
    eta_AFT = rep(eta_AFT, n),  # Same for all patients in this combo
    T_0 = T_0,
    points = points
  ))
}

apply_censoring <- function(points, censoring_time) {
  time_obs <- pmin(points, censoring_time)
  event_obs <- as.numeric(points <= censoring_time)
  censoring_rate <- 1 - mean(event_obs)
  
  return(list(
    time_obs = time_obs,
    event_obs = event_obs,
    censoring_rate = censoring_rate
  ))
}



###Run simulation
datasets <- list()
metadata_summary <- data.frame()

# Loop through sample sizes
for (n in n_points) {
  # Loop through all covariate combinations
  for (i in 1:nrow(covariate_combos)) {
    combo <- covariate_combos[i, ]
    combo_label <- combo$combo_label
    
    # Generate n_points event times for this covariate combination
    true_times <- generate_AFT_times(n, combo$age, combo$smoker, combo$sex,
                                     alpha, beta, b_age, b_smoker, b_sex)
    
    # Create base dataset with true times
    base_data <- data.frame(
      id = 1:n,
      age = rep(combo$age, n),
      smoker = rep(combo$smoker, n),
      sex = rep(combo$sex, n),
      combo_label = rep(combo_label, n),
      eta_AFT = true_times$eta_AFT,
      T_0 = true_times$T_0,
      points = true_times$points
    )
    
    # Apply each censoring scheme
    for (c_time in censoring_times) {
      c_name <- paste0("time", c_time)
      
      # Apply censoring
      censored <- apply_censoring(base_data$points, c_time)
      
      # Create censored dataset
      censored_data <- data.frame(
        base_data,
        censoring_time = c_time,
        time_obs = censored$time_obs,
        event_obs = censored$event_obs,
        censoring_rate = censored$censoring_rate
      )
      
      # Dataset name for flat structure
      dataset_name <- paste0("n", n, "_", combo_label, "_", c_name)
      
      # Store in flat structure
      datasets[[dataset_name]] <- list(
        data = censored_data,
        metadata = list(
          n_points = n,
          combo_label = combo_label,
          age = combo$age,
          smoker = combo$smoker,
          sex = combo$sex,
          censoring_time = c_time,
          censoring_rate = censored$censoring_rate,
          alpha = alpha,
          beta = beta,
          L_med = L_med,
          b_age = b_age,
          b_smoker = b_smoker,
          b_sex = b_sex,
          eta_AFT = true_times$eta_AFT[1],
          mean_points = mean(base_data$points),
          median_points = median(base_data$points),
          mean_time_obs = mean(censored$time_obs),
          n_events = sum(censored$event_obs),
          n_censored = sum(1 - censored$event_obs)
        )
      )
      
      # Add to summary table
      metadata_summary <- rbind(metadata_summary, data.frame(
        dataset_name = dataset_name,
        n_points = n,
        combo_label = combo_label,
        age = combo$age,
        smoker = combo$smoker,
        sex = combo$sex,
        censoring_time = c_time,
        censoring_rate = censored$censoring_rate,
        n_events = sum(censored$event_obs),
        n_censored = sum(1 - censored$event_obs),
        mean_points = mean(base_data$points),
        median_points = median(base_data$points),
        mean_time_obs = mean(censored$time_obs),
        median_time_obs = median(censored$time_obs),
        alpha = alpha,
        beta = beta,
        L_med = L_med,
        eta_AFT = true_times$eta_AFT[1]
      ))
    }
  }
}

#print num of datasets and metadata summary
length(datasets)
metadata_summary 
#contians summary stats about the datasets (could ignore - was for me to make sure all made sense)
#creates 72 datasets that includes all 8 covariate combinations, 3 censoring times and 3 sample sizes 
# for each dataset in datasets there are points (true values)
#time_obs (censored values)
#event_obs (is point was censored or not) 


#How to access datasets for analysis
#datasets[['dataset_name']] will give you the dataset
#censored points are found with datasets[['dataset_name']]$data$time_obs
#true points are found with datasets[['dataset_name']]$data$points
#find all datasets with names(datasets)



###Parameter Adjustment
library(dplyr)
library(ggplot2)

#median censoring_rate from different censoring_time and sample size
parameter_summary <- metadata_summary %>%
  group_by(n_points, censoring_time) %>%
  summarise(
    median_censoring_rate = median(censoring_rate),
    median_time_obs = median(median_time_obs),
    .groups = "drop"
  )
print(parameter_summary)

#ggplot for the parameters
ggplot(parameter_summary, aes(x = censoring_time, 
                              y = median_censoring_rate, 
                              color = factor(n_points))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Impact of Sample Size and Study Duration on Censoring Rate",
    subtitle = "Parameter Adjustment: varying sample size and censoring time",
    x = "Study Duration (years)",
    y = "Median Censoring Rate",
    color = "Sample Size"
  ) +
  theme_minimal(base_size = 14)

#median censoring_rate in all covariance combination
combo_summary <- metadata_summary %>%
  group_by(combo_label, n_points, censoring_time) %>%
  summarise(
    median_censoring_rate = median(censoring_rate),
    median_time_obs = median(median_time_obs),
    .groups = "drop"
  )

head(combo_summary)

#ggplot for all covariance combination
ggplot(combo_summary, aes(x = censoring_time, 
                          y = median_censoring_rate, 
                          color = combo_label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~n_points, ncol = 1) +
  labs(
    title = "Censoring Rate by Covariate Combination",
    x = "Study Duration (years)",
    y = "Median Censoring Rate",
    color = "Covariate Combination"
  ) +
  theme_minimal(base_size = 12)

###Model Comparison
library(survival)
library(purrr)

all_predictions <- list()

for (dname in names(datasets)) {
  dat <- datasets[[dname]]$data
  
  # Skip dataset if too few events
  if (sum(dat$event_obs) < 5) next
  
  # Train & Test set split
  set.seed(123)
  train_idx <- sample(1:nrow(dat), 0.7 * nrow(dat))
  train_data <- dat[train_idx, ]
  test_data <- dat[-train_idx, ]
  
  # Weibull AFT Model & Median Survival Time
  aft_fit <- survreg(Surv(time_obs, event_obs) ~ 1, 
                     data = train_data, dist = "weibull")
  aft_median_pred <- tryCatch(predict(aft_fit, newdata = test_data, type = "response"),
                              error = function(e) {rep(NA, nrow(test_data))})
  
  # Cox Model & Median Survival Time
  get_median_surv <- function(model, new_data_row) {fit <- survfit(model, newdata = new_data_row)
    median_time <- summary(fit)$table["median"]
    if (is.null(median_time) || is.na(median_time)) {
      return(NA)}
    return(unname(median_time))}
  
  cox_fit <- coxph(Surv(time_obs, event_obs) ~ 1, data = train_data)
  cox_median_pred <- tryCatch({purrr::map_dbl(1:nrow(test_data), 
                                              ~ get_median_surv(cox_fit, test_data[.x, , drop = FALSE]))}, error = function(e) {
                                                rep(NA, nrow(test_data))})
  
  # Naive Model & Median Survival Time
  naive_fit <- lm(time_obs ~ 1, data = train_data)
  
  naive_median_pred <- tryCatch(predict(naive_fit, newdata = test_data),
                                error = function(e) {rep(NA, nrow(test_data))})
  
  # All the predicted median value in here
  prediction_results <- test_data %>%
    dplyr::select(time_obs, event_obs, age, smoker, sex) %>%
    mutate(
      dataset = dname,
      n_points = datasets[[dname]]$metadata$n_points,
      censoring_time = datasets[[dname]]$metadata$censoring_time,
      aft_median_pred = aft_median_pred,
      cox_median_pred = cox_median_pred,
      naive_median_pred = naive_median_pred
    )
  
  all_predictions[[dname]] <- prediction_results
  print(all_predictions[[dname]])
}

# for further evaluation, the hate y is included in the [prediction_results]
# the real y for comparing in the test set named [test_data].
# hate y is predicted median values from these three model: AFT/COX/NATIVE
# some cox_median_pred is NA, due to the survival probability doesn't dip below 0.5 within the observed time frame of the Kaplan-Meier curve.
# thus, we can know that the Cox model does not cover all scenarios (returns NA). 
# in terms of predicting median time, while the AFT model can (always returns a number)- better in survival analysis.
