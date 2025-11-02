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
alpha <- 2.1                # Weibull shape parameter
L_med <- 8                  # Desired median lifetime for baseline group (years)
beta <- L_med / (log(2))^(1 / alpha)  # Weibull scale parameter

# AFT coefficients (for binary covariates)
b_age <- 0.02                 # Age effect 
b_smoker <- -0.5              # Smoking effect
b_sex <- 0.1                 # Sex effect

# All 8 covariate combinations, all binary covariates
generate_AFT_times_continuous <- function(n, alpha, beta, b_age, b_smoker, b_sex) {
  
  # Generate n continuous covariates for this dataset
  age_vals <- rnorm(n, mean = 70, sd = 10) # Age ~ N(70, 10)
  smoker_vals <- runif(n, 0, 1) # Smoking index ~ U(0, 1)
  sex_vals <- rbinom(n, 1, 0.5) # Binary, 0 = Female, 1 = Male (50% split)
  
  # Calculate η_AFT for *each* patient based on their unique covariates
  eta_AFT <- b_age * age_vals + b_smoker * smoker_vals + b_sex * sex_vals
  
  # Generate baseline event times
  T_0 <- rweibull(n, shape = alpha, scale = beta)
  
  # Calculate individual true event times
  points <- T_0 * exp(eta_AFT)
  
  return(list(
    age = age_vals,
    smoker = smoker_vals,
    sex = sex_vals,
    eta_AFT = eta_AFT,
    T_0 = T_0,
    points = points
  ))
}

#variable list
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
n_datasets <- 1000

# Loop through sample sizes
for (i in 1:n_datasets) {
  n <- sample(n_points, 1)
  c_time <- sample(censoring_times, 1)
  c_name <- paste0("time", c_time)
  
  #n_points event times with covariates
  true_times <- generate_AFT_times_continuous(n, alpha, beta, b_age, b_smoker, b_sex)
  
  #base dataset with true times and unique covariates
  base_data <- data.frame(
    id = 1:n,
    age = true_times$age,
    smoker = true_times$smoker,
    sex = true_times$sex,
    eta_AFT = true_times$eta_AFT,
    T_0 = true_times$T_0,
    points = true_times$points
  )
  
  #censoring
  censored <- apply_censoring(base_data$points, c_time)
  
  #censored dataset
  censored_data <- data.frame(
    base_data,
    censoring_time = c_time,
    time_obs = censored$time_obs,
    event_obs = censored$event_obs,
    censoring_rate = censored$censoring_rate
  )
  
  #dataset name
  dataset_name <- paste0("dataset_", i)
  
  #flat structure
  datasets[[dataset_name]] <- list(
    data = censored_data,
    metadata = list(
      dataset_name = dataset_name,
      n_points = n,
      censoring_time = c_time,
      censoring_rate = censored$censoring_rate,
      alpha = alpha,
      beta = beta,
      L_med = L_med,
      b_age = b_age,
      b_smoker = b_smoker,
      b_sex = b_sex,
      mean_age = mean(true_times$age),
      mean_smoker = mean(true_times$smoker),
      mean_sex = mean(true_times$sex),
      mean_eta_AFT = mean(true_times$eta_AFT),
      mean_points = mean(base_data$points),
      median_points = median(base_data$points),
      mean_time_obs = mean(censored$time_obs),
      median_time_obs = median(censored$time_obs),
      n_events = sum(censored$event_obs),
      n_censored = sum(1 - censored$event_obs)
    )
  )
  
  #summary table
  metadata_summary <- rbind(metadata_summary, data.frame(
    dataset_name = dataset_name,
    n_points = n,
    censoring_time = c_time,
    censoring_rate = censored$censoring_rate,
    n_events = sum(censored$event_obs),
    n_censored = sum(1 - censored$event_obs),
    mean_age = mean(true_times$age),
    mean_smoker = mean(true_times$smoker),
    mean_sex = mean(true_times$sex),
    mean_eta_AFT = mean(true_times$eta_AFT),
    mean_points = mean(base_data$points),
    median_points = median(base_data$points),
    mean_time_obs = mean(censored$time_obs),
    median_time_obs = median(censored$time_obs),
    alpha = alpha,
    beta = beta,
    L_med = L_med
    # Note: 'eta_AFT' is no longer a single value, so we store its mean
  ))
}


#print num of datasets and metadata summary
cat("Total datasets generated:", length(datasets), "\n")
print("Head of metadata summary:")
head(metadata_summary)

#contians summary stats about the datasets (could ignore - was for me to make sure all made sense)
#creates 1000 datasets that includes all 8 covariate combinations
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
  aft_fit <- survreg(Surv(time_obs, event_obs) ~ age + smoker + sex, 
                     data = train_data, dist = "weibull")
  aft_median_pred <- tryCatch(predict(aft_fit, newdata = test_data, type = "response"),
                              error = function(e) {rep(NA, nrow(test_data))})
  
  # Cox Model & Median Survival Time
  get_median_surv <- function(model, new_data_row) {fit <- survfit(model, newdata = new_data_row)
                      median_time <- summary(fit)$table["median"]
                      if (is.null(median_time) || is.na(median_time)) {
                        return(NA)}
                      return(unname(median_time))}
  
  cox_fit <- coxph(Surv(time_obs, event_obs) ~ age + smoker + sex, data = train_data)
  cox_median_pred <- tryCatch({purrr::map_dbl(1:nrow(test_data), 
                                              ~ get_median_surv(cox_fit, test_data[.x, , drop = FALSE]))}, error = function(e) {
                                                rep(NA, nrow(test_data))})
  
  # Naive Model & Median Survival Time
  naive_fit <- lm(time_obs ~ age + smoker + sex, data = train_data)
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
