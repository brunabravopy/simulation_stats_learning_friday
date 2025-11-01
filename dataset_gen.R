# AFT Weibull Simulation with Multiple Sample Sizes and Censoring Schemes
# Model: T = T_0 * exp(η_AFT) where T_0 ~ Weibull(shape=α, scale=β)
# η_AFT = b_age * age_binary + b_smoker * smoker + b_sex * sex
# Covariates: age (0=<70, 1=≥70), smoker (0=non-smoker, 1=smoker), sex (0=female, 1=male)
#
# For each covariate combination, generate n_points patients from that distribution
# Then apply 3 censoring times to each dataset

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

#run simulation
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
metadata_summary #contians summary stats about the datasets (could ignore - was for me to make sure all made sense)
#creates 72 datasets that includes all 8 covariate combinations, 3 censoring times and 3 sample sizes 
# for each dataset in datasets there are points (true values)
    #time_obs (censored values)
    #event_obs (is point was censored or not) 


#How to access datasets for analysis
    #datasets[['dataset_name']] will give you the dataset
    #censored points are found with datasets[['dataset_name']]$data$time_obs
    #true points are found with datasets[['dataset_name']]$data$points
#find all datasets with names(datasets)