# Life-time simulation (Weibull + Cox-like effects)

# Packages
set.seed(1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)})

# Output folder on Desktop 
get_desktop <- function() {
  # Windows
  dp_win <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
  if (nzchar(Sys.getenv("USERPROFILE")) && dir.exists(dp_win)) return(dp_win)
  # Fallback: home
  return(path.expand("~"))}

desktop_path <- get_desktop()
out_dir <- file.path(desktop_path, "Weibull_Survival_Plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Saving plots to:\n", out_dir, "\n\n", sep = "")

# Base parameters
n <- 1000
alpha <- 3.5                # shape: risk grows with age
start_age_mean <- 30        # average age at start
L_med <- 45                 # desired median of life time (years after start)

# Auto-scale (beta) to hit chosen median
beta <- L_med / (log(2))^(1 / alpha)

# individual variables
start_age <- pmin(pmax(rnorm(n, start_age_mean, 5), 20), 40)
sex <- sample(0:1, n, replace = TRUE)                   # 0=female, 1=male
fitness <- pmin(pmax(rnorm(n, 0.5, 0.15), 0), 1)        # 0–1
effectiveness <- pmin(pmax(rnorm(n, 0.6, 0.2), 0), 1)   # 0–1

# Cox-like coefficients 
coef_age <- log(1.2) / 10       # +20% risk per +10y start age
coef_sex <- log(1.05)           # males +5% risk
coef_fit <- log(0.97) / 0.1     # +0.1 fitness -> ~-3% risk
coef_eff <- log(0.7)            # effectiveness -> ~-30% risk

risk_score <- coef_age * start_age + coef_sex * sex + coef_fit * fitness + coef_eff * effectiveness

# Generate life time Weibull
U <- runif(n)
time_to_death <- beta * ((-log(U)) / exp(risk_score))^(1 / alpha)
final_age <- start_age + time_to_death


# Small KM helper (no censoring here)
calc_surv <- function(time_vec) {
  fit <- survfit(Surv(time_vec, rep(1, length(time_vec))) ~ 1)
  data.frame(time = fit$time, S = fit$surv, stringsAsFactors = FALSE)
}

# Build main data frame
df <- data.frame(
  start_age = start_age,
  final_age = final_age,
  time = time_to_death,
  sex = factor(sex, labels = c("Female", "Male")),
  fitness = fitness,
  effectiveness = effectiveness
)


# KM by SEX
km_female <- calc_surv(df$final_age[df$sex == "Female"])
km_male   <- calc_surv(df$final_age[df$sex == "Male"])

p1 <- ggplot() +
  geom_step(data = km_female, aes(x = time, y = S, color = "Female")) +
  geom_step(data = km_male,   aes(x = time, y = S, color = "Male")) +
  scale_color_manual(values = c("Female" = "magenta", "Male" = "blue")) +
  labs(x = "Age (years)", y = "Alive proportion (survival)",
       title = "Kaplan–Meier by sex (alpha = 3.5)", color = "Group") +
  theme_minimal()

p2 <- ggplot() +
  geom_step(data = km_female, aes(x = time, y = 1 - S, color = "Female")) +
  geom_step(data = km_male,   aes(x = time, y = 1 - S, color = "Male")) +
  scale_color_manual(values = c("Female" = "magenta", "Male" = "blue")) +
  labs(x = "Age (years)", y = "Cumulative probability of death",
       title = "Cumulative death by sex (alpha = 3.5)", color = "Group") +
  theme_minimal()

ggsave(file.path(out_dir, "KM_by_sex_R.png"), p1, dpi = 300, width = 8, height = 5)
ggsave(file.path(out_dir, "CDF_by_sex_R.png"), p2, dpi = 300, width = 8, height = 5)

# KM by FITNESS level
fit_cat <- cut(df$fitness, breaks = c(0, 0.4, 0.6, 1),
               labels = c("Low", "Medium", "High"), include.lowest = TRUE)
df$fit_cat <- fit_cat

# Build KM table per category 
split_fit <- split(df$final_age, df$fit_cat)
km_fit_list <- lapply(names(split_fit), function(level) {
  tmp <- calc_surv(split_fit[[level]])
  tmp$fit_cat <- level
  tmp
})
km_fit <- do.call(rbind, km_fit_list)

p3 <- ggplot(km_fit, aes(x = time, y = S, color = fit_cat)) +
  geom_step() +
  scale_color_manual(values = c("Low" = "red", "Medium" = "orange", "High" = "green")) +
  labs(x = "Age (years)", y = "Alive proportion (survival)",
       title = "Kaplan–Meier by fitness (alpha = 3.5)", color = "Fitness") +
  theme_minimal()

p4 <- ggplot(km_fit, aes(x = time, y = 1 - S, color = fit_cat)) +
  geom_step() +
  scale_color_manual(values = c("Low" = "red", "Medium" = "orange", "High" = "green")) +
  labs(x = "Age (years)", y = "Cumulative probability of death",
       title = "Cumulative death by fitness (alpha = 3.5)", color = "Fitness") +
  theme_minimal()

ggsave(file.path(out_dir, "KM_by_fitness_R.png"), p3, dpi = 300, width = 8, height = 5)
ggsave(file.path(out_dir, "CDF_by_fitness_R.png"), p4, dpi = 300, width = 8, height = 5)


# 9) KM by EFFECTIVENESS level
eff_cat <- cut(df$effectiveness, breaks = c(0, 0.4, 0.7, 1),
               labels = c("Low", "Medium", "High"), include.lowest = TRUE)
df$eff_cat <- eff_cat

split_eff <- split(df$final_age, df$eff_cat)
km_eff_list <- lapply(names(split_eff), function(level) {
  tmp <- calc_surv(split_eff[[level]])
  tmp$eff_cat <- level
  tmp
})
km_eff <- do.call(rbind, km_eff_list)

p5 <- ggplot(km_eff, aes(x = time, y = S, color = eff_cat)) +
  geom_step() +
  scale_color_manual(values = c("Low" = "red", "Medium" = "orange", "High" = "green")) +
  labs(x = "Age (years)", y = "Alive proportion (survival)",
       title = "Kaplan–Meier by treatment effectiveness (alpha = 3.5)",
       color = "Effectiveness") +
  theme_minimal()

p6 <- ggplot(km_eff, aes(x = time, y = 1 - S, color = eff_cat)) +
  geom_step() +
  scale_color_manual(values = c("Low" = "red", "Medium" = "orange", "High" = "green")) +
  labs(x = "Age (years)", y = "Cumulative probability of death",
       title = "Cumulative death by effectiveness (alpha = 3.5)",
       color = "Effectiveness") +
  theme_minimal()

ggsave(file.path(out_dir, "KM_by_effectiveness_R.png"), p5, dpi = 300, width = 8, height = 5)
ggsave(file.path(out_dir, "CDF_by_effectiveness_R.png"), p6, dpi = 300, width = 8, height = 5)


#Histogram + density of final ages
# -------------------------
p7 <- ggplot(df, aes(x = final_age)) +
  geom_histogram(aes(y = ..density..), bins = 40,
                 fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1.2) +
  labs(x = "Age at death (years)", y = "Probability density",
       title = "Simulated distribution of ages at death (Weibull alpha = 3.5)") +
  theme_minimal()

ggsave(file.path(out_dir, "Histogram_final_age_alpha3.5_R.png"),
       p7, dpi = 300, width = 9, height = 5)

#  Final summary
png_files <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
cat("\nSaved PNG files:\n", paste0(" - ", basename(png_files), collapse = "\n"), "\n", sep = "")
cat(sprintf("\nWeibull scale (beta): %.2f\n", beta))
cat(sprintf("Mean start age: %.1f\n", mean(start_age)))
cat(sprintf("Median age at death (Females): %.1f\n", median(final_age[sex == 0])))
cat(sprintf("Median age at death (Males): %.1f\n", median(final_age[sex == 1])))
cat(sprintf("Max age observed: %.1f\n", max(final_age)))
