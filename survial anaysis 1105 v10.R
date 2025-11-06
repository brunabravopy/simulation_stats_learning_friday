# AFT Weibull Simulation with Multiple Sample Sizes and Censoring Schemes
# Model: T = T_0 * exp(η_AFT) where T_0 ~ Weibull(shape=α, scale=β)
# η_AFT = b_age * age_binary + b_smoker * smoker + b_sex * sex

library(dplyr)
library(tidyr)
library(survival)
library(purrr)
library(ggplot2)
library(patchwork)

### 1. Generate the data
set.seed(3456)

# Sample sizes to generate
n_points <- c(100, 300, 600)

# Censoring times (study end dates) #Chronic Non-Communicable Diseases
time1 <- 4 # years
time2 <- 6 # years  
time3 <- 10 # years
censoring_times <- c(time1, time2, time3)

# AFT Weibull parameters
alpha <- 3 # Weibull shape parameter
L_med <- 6  # Desired median lifetime for baseline group (years)
beta <- L_med / (log(2))^(1 / alpha)  # Weibull scale parameter

# AFT coefficients
b_age <- -0.25                 # Age effect
b_age2   <- 0.15 
b_smoker <- -0.6              # Smoking effect
b_sex <- 0.2                 # Sex effect

# Generate the data
gen_data <- function(n) {
  age_raw <- rnorm(n, 70, 10) #mean = 70, sd = 10
  age_std <- (age_raw - 70) / 10
  smoker  <- rbinom(n, 1, 0.3) #binary, with 30% split randomly (30% = 1 = smoker)
  sex     <- rbinom(n, 1, 0.5) #binary, with 50% split randomly
  
  # Model with parameters  
  lp <- b_age * age_std + b_age2 * (age_std^2) + b_smoker * smoker + b_sex * sex
  
  # Baseline Weibull
  T0 <- rweibull(n, shape = alpha, scale = beta)
  
  # Model with weibull baseline
  true_time <- T0 * exp(lp)
  
  # Data frame
  data.frame(
    age_std  = age_std,
    age_std2 = age_std^2,
    smoker   = smoker,
    sex      = sex,
    true_time = true_time
  )
}

# Apply censoring -> ctime = study time
apply_censoring <- function(true_time, ctime) {
  time_obs  <- pmin(true_time, ctime)
  event_obs <- as.numeric(true_time <= ctime)
  data.frame(
    time_obs = time_obs,
    event_obs = event_obs,
    censoring_rate = 1 - mean(event_obs) # mean() -> percent of event happening
  )
}

# Show an example of generated data
set.seed(123)
#example - using actual simulation parameters
n <- 300  # one of the sample sizes used: c(100, 300, 600)
ctime <- 6  # one of the censoring times used: c(4, 6, 10)
#from generated data
example_df <- gen_data(n)
censored_df <- apply_censoring(example_df$true_time, ctime)
example_data <- cbind(example_df, censored_df)

# Table for example
print(head(example_data, 10))

# p1: Age vs True survival time
p1 <- ggplot(example_data, aes(x = age_std, y = true_time, color = factor(smoker))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("forestgreen", "red3"),
                     labels = c("Non-smoker", "Smoker")) +
  labs(
    title = "Age and Smoking Effects on True Survival",
    x = "Standardized Age", y = "True Survival Time (years)",
    color = "Smoking Status"
  ) +
  theme_minimal(base_size = 13)

# p2: True vs Observed survival times
p2 <- ggplot(example_data, aes(x = true_time)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", alpha = 0.6) +
  geom_histogram(aes(x = time_obs), binwidth = 0.5, fill = "tomato", alpha = 0.5) +
  labs(
    title = "True vs Observed Survival Times",
    x = "Time (years)", y = "Count",
    caption = "Blue = True times; Red = Observed (after censoring)"
  ) +
  theme_minimal(base_size = 13)

# Simulation Overview
overview_plot <- (p1 | p2) +
  plot_annotation(
    title = "Simulation Overview",
    subtitle = "Example data from AFT-Weibull model with censoring (study time = 6 years / sample size = 300).",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

overview_plot



### 2. Run simulation
n_datasets <- 100
datasets <- list()

for (i in 1:n_datasets) {
  n     <- sample(n_points, 1)
  ctime <- sample(censoring_times, 1)
  df  <- gen_data(n)
  cen <- apply_censoring(df$true_time, ctime)
  datasets[[paste0("dataset_", i)]] <- cbind(
    df,
    cen,
    n_points = n,
    censoring_time = ctime
  )
}

cat("datasets:", length(datasets), "\n")



### 3. Model Comparison
all_preds <- list()

for (dname in names(datasets)) {
  dat <- datasets[[dname]]
  if (sum(dat$event_obs) < 5) next
  
  # Split trainning set and test set
  set.seed(123)
  idx   <- sample(1:nrow(dat), floor(0.7 * nrow(dat))) # 0.7/0.3
  train <- dat[idx, ]
  test  <- dat[-idx, ]
  
  ## AFT model (Weibull)
  aft_fit <- tryCatch(
    survreg(
      Surv(time_obs, event_obs) ~ age_std + smoker + sex,
      data = train,
      dist = "weibull"
    ),
    error = function(e) NULL
  )
  
  aft_pred <- if (!is.null(aft_fit)) {
    predict(aft_fit, newdata = test, type = "quantile", p = 0.5)
  } else rep(NA_real_, nrow(test)) # predict the median survival time
  
  aft_pred_cap <- pmin(aft_pred, 25)  # tighter numerical cap
  
  ## Cox model
  cox_fit <- tryCatch(
    coxph(Surv(time_obs, event_obs) ~ age_std + smoker + sex, data = train),
    error = function(e) NULL
  )
  
  get_median_from_cox <- function(model, newrow) {
    bh <- basehaz(model, centered = FALSE)
    lp <- as.numeric(predict(model, newdata = newrow, type = "lp"))
    s  <- exp(-bh$hazard * exp(lp))
    idx <- which(s <= 0.5)[1]
    if (is.na(idx)) return(NA_real_)
    bh$time[idx]
  }
  
  cox_pred <- if (!is.null(cox_fit)) {
    purrr::map_dbl(1:nrow(test), ~ get_median_from_cox(cox_fit, test[.x, , drop = FALSE]))
  } else rep(NA_real_, nrow(test)) # predict the median survival time
  
  cox_pred_cap <- pmin(cox_pred, 25) # tighter numerical cap
  
  ## Naive model
  naive_fit <- lm(time_obs ~ age_std + smoker + sex, data = train)
  naive_pred <- predict(naive_fit, newdata = test)
  naive_pred_cap <- pmin(naive_pred, 25)
  
  ## All predictor
  all_preds[[dname]] <- data.frame(
    dataset_name = dname,
    true_time = test$true_time,
    aft_raw   = aft_pred_cap,
    cox_raw   = cox_pred_cap,
    naive_raw = naive_pred_cap
  )
}

# Show predictors
res_df <- bind_rows(all_preds, .id = "dataset_name")

res_long_raw <- res_df %>%
  select(dataset_name, true_time, aft_raw, cox_raw, naive_raw) %>%
  pivot_longer(cols = starts_with(c("aft", "cox", "naive")),
               names_to = "model", values_to = "pred_time") %>%
  mutate(model = factor(model,
                        levels = c("aft_raw", "cox_raw", "naive_raw"),
                        labels = c("AFT (Weibull)", "Cox", "Naive")))

# Predicted vs True Scatterplot
ggplot(res_long_raw, aes(x = true_time, y = pred_time, color = model)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~model, scales = "free", nrow = 1) +
  labs(
    title = "Predicted vs True Survival Times",
    subtitle = "Model comparison across simulated datasets",
    x = "True Survival Time (years)",
    y = "Predicted Median Survival Time (years)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 25))

# Predictor with different parameters
res_aft <- bind_rows(all_preds, .id = "dataset_name") %>%
  select(dataset_name, true_time, aft_pred = aft_raw) %>%
  left_join(
    bind_rows(datasets, .id = "dataset_name") %>%
      select(dataset_name, n_points, censoring_time),
    by = "dataset_name"
  ) %>%
  mutate(
    abs_error = abs(aft_pred - true_time),
    bias = aft_pred - true_time
  )

summary_aft <- res_aft %>%
  group_by(n_points, censoring_time) %>%
  summarise(
    mean_pred = mean(aft_pred, na.rm = TRUE),
    mean_true = mean(true_time, na.rm = TRUE),
    .groups = "drop"
  )

summary_aft_long <- summary_aft %>%
  pivot_longer(cols = c(mean_pred, mean_true),
               names_to = "Type", values_to = "Mean_Survival") %>%
  mutate(
    Type = recode(Type,
                  "mean_pred" = "Predicted (AFT)",
                  "mean_true" = "True (Simulated)")
  )


### 4. Evaluation (Median-based approach)
res <- bind_rows(all_preds)

##Winsorize helper: trim 1% extremes to avoid outliers
winsorize_vec <- function(x, lower_p = 0.01, upper_p = 0.99) {
  q <- quantile(x, probs = c(lower_p, upper_p), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

##Prepare long format for median-based computation
res_long_raw <- res %>%
  pivot_longer(cols = c(aft_raw, cox_raw, naive_raw),
               names_to = "model", values_to = "pred")

##For each dataset & model: compute median(pred) & median(true_time)
median_summary <- res_long_raw %>%
  group_by(dataset_name, model) %>%
  summarise(
    median_pred = median(winsorize_vec(pred), na.rm = TRUE),
    median_true = median(true_time, na.rm = TRUE),  # True values not winsorized
    .groups = "drop"
  )

##Add censoring_time and n_points information
median_summary <- median_summary %>%
  left_join(
    bind_rows(datasets, .id = "dataset_name") %>%
      select(dataset_name, censoring_time, n_points) %>%
      distinct(),
    by = "dataset_name"
  )

##Compute squared error and bias between medians
median_summary <- median_summary %>%
  mutate(
    mse_median = (median_pred - median_true)^2,
    bias_median = median_pred - median_true
  )

##Aggregate across datasets (main metrics) - overall and by study length
overall_median <- median_summary %>%
  group_by(model) %>%
  summarise(
    Median_MSE = mean(mse_median, na.rm = TRUE),
    Median_Bias = mean(bias_median, na.rm = TRUE),
    .groups = "drop"
  )

##Metrics stratified by study length
median_by_study <- median_summary %>%
  group_by(model, censoring_time) %>%
  summarise(
    Median_MSE = mean(mse_median, na.rm = TRUE),
    Median_Bias = mean(bias_median, na.rm = TRUE),
    .groups = "drop"
  )

##Metrics stratified by sample size
median_by_n <- median_summary %>%
  group_by(model, n_points) %>%
  summarise(
    Median_MSE = mean(mse_median, na.rm = TRUE),
    Median_Bias = mean(bias_median, na.rm = TRUE),
    .groups = "drop"
  )

cat("MEDIAN-BASED METRICS\n")
print(overall_median)
cat("\nMEDIAN-BASED METRICS BY STUDY LENGTH\n")
print(median_by_study)
cat("\nMEDIAN-BASED METRICS BY SAMPLE SIZE\n")
print(median_by_n)


##CI Coverage (Median-based)
# For each model, censoring_time, AND n_points, compute CI range based on median_pred
# This ensures CI is computed separately for each combination of censoring scenario and sample size
ci_ranges <- median_summary %>%
  group_by(model, censoring_time, n_points) %>%
  summarise(
    ci_low = quantile(median_pred, 0.025, na.rm = TRUE),
    ci_high = quantile(median_pred, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Check if median_true falls within CI for each dataset
coverage_summary <- median_summary %>%
  left_join(ci_ranges, by = c("model", "censoring_time", "n_points")) %>%
  mutate(
    covered = as.numeric(median_true >= ci_low & median_true <= ci_high)
  )

# Calculate mean coverage across datasets for each model (overall and by study length)
overall_coverage <- coverage_summary %>%
  group_by(model) %>%
  summarise(
    Mean_Coverage_95 = mean(covered, na.rm = TRUE),
    .groups = "drop"
  )

##Coverage stratified by study length
coverage_by_study <- coverage_summary %>%
  group_by(model, censoring_time) %>%
  summarise(
    Mean_Coverage_95 = mean(covered, na.rm = TRUE),
    .groups = "drop"
  )

##Coverage stratified by sample size
coverage_by_n <- coverage_summary %>%
  group_by(model, n_points) %>%
  summarise(
    Mean_Coverage_95 = mean(covered, na.rm = TRUE),
    .groups = "drop"
  )

cat("95% CI COVERAGE (Median-based)\n")
print(overall_coverage)
cat("\n95% CI COVERAGE BY STUDY LENGTH\n")
print(coverage_by_study)
cat("\n95% CI COVERAGE BY SAMPLE SIZE\n")
print(coverage_by_n)


### 5. Visualization of Model Performance (Updated for Median-based Evaluation)
## Rename models for clarity
overall_median$model <- recode(overall_median$model,
                               "aft_raw" = "AFT (Weibull)",
                               "cox_raw" = "Cox",
                               "naive_raw" = "Naive")

overall_coverage$model <- recode(overall_coverage$model,
                                 "aft_raw" = "AFT (Weibull)",
                                 "cox_raw" = "Cox",
                                 "naive_raw" = "Naive")

## Combine all metrics into one unified table
measure_summary <- overall_median %>%
  rename(MSE = Median_MSE, Bias = Median_Bias) %>%
  left_join(overall_coverage, by = "model")

cat("\n==== Combined Model Performance Table (Median-based Evaluation) ====\n")
print(measure_summary)

## Reshape for plotting
measure_long <- measure_summary %>%
  pivot_longer(cols = -model, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric_Group = case_when(
      grepl("MSE", Metric) ~ "Median-based Mean Squared Error",
      grepl("Bias", Metric) ~ "Median-based Bias",
      grepl("Coverage", Metric) ~ "95% CI Coverage (Median-based)"
    )
  )

## Colors and plot
metric_colors <- c(
  "MSE" = "#4E79A7",
  "Bias" = "#59A14F",
  "Mean_Coverage_95" = "#76B7B2"
)


## Stratified visualization by study length
# Rename models for stratified plot
median_by_study$model <- recode(median_by_study$model,
                               "aft_raw" = "AFT (Weibull)",
                               "cox_raw" = "Cox",
                               "naive_raw" = "Naive")

coverage_by_study$model <- recode(coverage_by_study$model,
                                 "aft_raw" = "AFT (Weibull)",
                                 "cox_raw" = "Cox",
                                 "naive_raw" = "Naive")

# Combine MSE, Bias, and Coverage for stratified plot
stratified_data <- median_by_study %>%
  pivot_longer(cols = c(Median_MSE, Median_Bias),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = recode(Metric,
                   "Median_MSE" = "MSE",
                   "Median_Bias" = "Bias"),
    censoring_time = factor(censoring_time)
  ) %>%
  bind_rows(
    coverage_by_study %>%
      rename(Value = Mean_Coverage_95) %>%
      mutate(Metric = "CI Coverage",
             censoring_time = factor(censoring_time))
  )

stratified_plot <- ggplot(stratified_data, 
                          aes(x = censoring_time, y = Value, 
                              color = model, group = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  geom_hline(
    data = subset(stratified_data, Metric == "Bias"),
    aes(yintercept = 0),
    color = "black", linetype = "dashed"
  ) +
  geom_hline(
    data = subset(stratified_data, Metric == "CI Coverage"),
    aes(yintercept = 0.95),
    color = "red", linetype = "dashed", linewidth = 1
  ) +
  labs(
    title = "Model Performance by Study Length (Censoring Time)",
    subtitle = "How model performance changes with different study durations",
    x = "Study Length / Censoring Time (years)",
    y = "Value",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

stratified_plot

## Stratified visualization by sample size
# Rename models for stratified plot by sample size
median_by_n$model <- recode(median_by_n$model,
                               "aft_raw" = "AFT (Weibull)",
                               "cox_raw" = "Cox",
                               "naive_raw" = "Naive")

coverage_by_n$model <- recode(coverage_by_n$model,
                                 "aft_raw" = "AFT (Weibull)",
                                 "cox_raw" = "Cox",
                                 "naive_raw" = "Naive")

# Combine MSE, Bias, and Coverage for stratified plot by sample size
stratified_data_n <- median_by_n %>%
  pivot_longer(cols = c(Median_MSE, Median_Bias),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = recode(Metric,
                   "Median_MSE" = "MSE",
                   "Median_Bias" = "Bias"),
    n_points = factor(n_points)
  ) %>%
  bind_rows(
    coverage_by_n %>%
      rename(Value = Mean_Coverage_95) %>%
      mutate(Metric = "CI Coverage",
             n_points = factor(n_points))
  )

stratified_plot_n <- ggplot(stratified_data_n, 
                          aes(x = n_points, y = Value, 
                              color = model, group = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  geom_hline(
    data = subset(stratified_data_n, Metric == "Bias"),
    aes(yintercept = 0),
    color = "black", linetype = "dashed"
  ) +
  geom_hline(
    data = subset(stratified_data_n, Metric == "CI Coverage"),
    aes(yintercept = 0.95),
    color = "red", linetype = "dashed", linewidth = 1
  ) +
  labs(
    title = "Model Performance by Sample Size",
    subtitle = "How model performance changes with different sample sizes",
    x = "Sample Size (n)",
    y = "Value",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

stratified_plot_n
