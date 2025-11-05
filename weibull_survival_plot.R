# Weibull Survival Time Plot Generator
# Helper script to visualize Weibull survival function with specified parameters

library(ggplot2)

# AFT Weibull parameters (from survial anaysis 1103 v4.R)
alpha <- 3 # Weibull shape parameter
L_med <- 6  # Desired median lifetime for baseline group (years)
beta <- L_med / (log(2))^(1 / alpha)  # Weibull scale parameter

# Generate time points for plotting
time_points <- seq(0, 20, by = 0.1)

# Calculate Weibull survival function: S(t) = exp(-(t/beta)^alpha)
survival_prob <- exp(-(time_points / beta)^alpha)

# Create data frame for plotting
plot_data <- data.frame(
  time = time_points,
  survival = survival_prob
)

# Plot: Survival Function
survival_plot <- ggplot(plot_data, aes(x = time, y = survival)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_vline(xintercept = L_med, linetype = "dashed", color = "red", linewidth = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1.5) +
  annotate("text", x = L_med + 1, y = 0.55, label = paste("Median =", L_med, "years"), 
           color = "red", size = 5, fontface = "bold") +
  labs(
    title = "Weibull Survival Function",
    subtitle = paste("Shape (α) =", alpha, ", Scale (β) =", round(beta, 3)),
    x = "Time (years)",
    y = "Survival Probability S(t)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5))


print(survival_plot)
