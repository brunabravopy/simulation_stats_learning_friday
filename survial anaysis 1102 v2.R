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
time1 <- 6 # years
time2 <- 8 # years  
time3 <- 10 # years
censoring_times <- c(time1, time2, time3)

# AFT Weibull parameters
alpha <- 3 # Weibull shape parameter
L_med <- 6  # Desired median lifetime for baseline group (years)
beta <- L_med / (log(2))^(1 / alpha)  # Weibull scale parameter

# AFT coefficients (for binary covariates)
b_age <- 0.25                 # Age effect
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
  lp <- b_age1 * age_std + b_age2 * (age_std^2) + b_smoker * smoker + b_sex * sex
    
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

# show an example of generated data
set.seed(123)
  #example
  n <- 300
  ctime <- 8
  #from generated data
  example_df <- gen_data(n)
  censored_df <- apply_censoring(example_df$true_time, ctime)
  example_data <- cbind(example_df, censored_df)

# table for example
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
    subtitle = "Data generated from AFT-Weibull model with censoring (study time = 8 years / sample size = 300)",
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
print("Head of metadata summary:")
head(metadata_summary)



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
      Surv(time_obs, event_obs) ~ age_std + I(age_std^2) + smoker + sex,
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
  
  ## all predictor
  all_preds[[dname]] <- data.frame(
    dataset_name = dname,
    true_time = test$true_time,
    aft_raw   = aft_pred_cap,
    cox_raw   = cox_pred_cap,
    naive_raw = naive_pred_cap,
    aft_log   = log(pmax(aft_pred, 1e-6)),
    cox_log   = log(pmax(cox_pred, 1e-6)),
    naive_log = log(pmax(naive_pred, 1e-6)),
    true_log  = log(test$true_time)
  )
}

# show predictors
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


### 4. Evaluation
res <- bind_rows(all_preds)

## small helper: winsorize to ignore 1% extreme values
winsorize_vec <- function(x, lower_p = 0.01, upper_p = 0.99) {
  q <- quantile(x, probs = c(lower_p, upper_p), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

## raw-scale MSE
res_long_raw <- res %>%
  pivot_longer(cols = c(aft_raw, cox_raw, naive_raw),
               names_to = "model", values_to = "pred")

raw_summary <- res_long_raw %>%
  group_by(dataset_name, model) %>%
  summarise(
    mse = mean( (winsorize_vec(pred) - winsorize_vec(true_time))^2, na.rm = TRUE ),
    med_diff = median(pred - true_time, na.rm = TRUE),
    .groups = "drop"
  )

overall_raw <- raw_summary %>%
  group_by(model) %>%
  summarise(
    Mean_MSE = mean(mse, na.rm = TRUE),
    Mean_MedDiff = mean(med_diff, na.rm = TRUE)
  )

# row mse table
cat("==== RAW-SCALE METRICS ====\n")
print(overall_raw)

## log-scale MSE
res_long_log <- res %>%
  pivot_longer(cols = c(aft_log, cox_log, naive_log),
               names_to = "model", values_to = "pred_log")

log_summary <- res_long_log %>%
  group_by(dataset_name, model) %>%
  summarise(
    mse_log = mean((pred_log - true_log)^2, na.rm = TRUE),
    med_diff_log = median(pred_log - true_log, na.rm = TRUE),
    .groups = "drop"
  )

overall_log <- log_summary %>%
  group_by(model) %>%
  summarise(
    Mean_MSE_log = mean(mse_log, na.rm = TRUE),
    Mean_MedDiff_log = mean(med_diff_log, na.rm = TRUE)
  )

# log mse table
cat("==== LOG-SCALE METRICS ====\n")
print(overall_log)

