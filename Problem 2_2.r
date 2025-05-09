# Load libraries
library(tidyverse)
library(lubridate)
library(patchwork)

# Read and prepare the data
setwd("C:/Users/snehi/Documents/Time Series Analysis/Assignment4")
data <- read_csv("transformer_data.csv")

# Rename columns for plotting
data <- data %>%
  rename(
    Y_t     = Y,
    T_a_t   = Ta,
    Phi_s_t = S,
    Phi_l_t = I
  )

# Plot all series
p1 <- ggplot(data, aes(x = time, y = Y_t)) +
  geom_line(color = "red") + labs(title = "Transformer Temperature", y = "Y_t (°C)", x = "")

p2 <- ggplot(data, aes(x = time, y = T_a_t)) +
  geom_line(color = "blue") + labs(title = "Outdoor Air Temperature", y = "T_a_t (°C)", x = "")

p3 <- ggplot(data, aes(x = time, y = Phi_s_t)) +
  geom_line(color = "orange") + labs(title = "Solar Radiation", y = "Phi_s_t (W/m²)", x = "")

p4 <- ggplot(data, aes(x = time, y = Phi_l_t)) +
  geom_line(color = "green") + labs(title = "Transformer Load", y = "Phi_l_t (kA)", x = "Time Step")

# Display combined plot
(p1 / p2 / p3 / p4) + plot_layout(guides = "collect")

# Prepare data for Kalman filter
df <- data %>%
  transmute(
    Y  = Y_t,
    Ta = T_a_t,
    S  = Phi_s_t,
    I  = Phi_l_t
  )

# === Kalman filter log-likelihood function ===
run_kalman_filter <- function(par, df) {
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  Sigma1lt <- matrix(exp(par[5]), 1, 1)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  C   <- matrix(1, 1, 1)
  Sigma2 <- matrix(exp(par[6]), 1, 1)
  X0  <- matrix(par[7], 1, 1)

  obs_cols <- c("Y")
  input_cols <- c("Ta", "S", "I")

  Y  <- as.matrix(df[, obs_cols])
  U  <- as.matrix(df[, input_cols])
  Tn <- nrow(df)

  n      <- nrow(A)
  x_est  <- X0
  P_est  <- diag(1e1, n)

  predictions <- numeric(Tn)
  residuals   <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t,], ncol=1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred <- C %*% x_pred
    S_t    <- C %*% P_pred %*% t(C) + Sigma2
    innov  <- Y[t, ] - y_pred

    # Save outputs
    predictions[t] <- y_pred
    residuals[t]   <- innov

    K_t   <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
  }

  list(predictions = predictions, residuals = residuals)
}


# === Estimation wrapper ===
estimate_dt <- function(start_par, df, lower=NULL, upper=NULL) {
  negLL <- function(par){ -kf_logLik_dt(par, df) }
  optim(
    par     = start_par,
    fn      = negLL,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit=1000, trace=1)
  )
}

# === Run estimation ===
start_par <- c(0.5, 0.1, 0.1, 0.1, log(0.5), log(0.5), 20)
lower     <- c(-1, -1, -1, -1, log(1e-6), log(1e-6), -100)
upper     <- c(1, 1, 1, 1, log(100), log(100), 100)

result <- estimate_dt(start_par, df, lower, upper)

# === Display results ===
print(result$par)

# Generating Plots
library(ggplot2)
library(gridExtra)
library(forecast)

# Run Kalman filter with estimated parameters
kf_output <- run_kalman_filter(result$par, df)

# Residuals
residuals <- kf_output$residuals

# Plot residual time series
p_res <- ggplot(data.frame(t = 1:length(residuals), resid = residuals), aes(x = t, y = resid)) +
  geom_line(color = "darkred") + theme_minimal() +
  labs(title = "Residuals", y = "Residual", x = "Time")

# ACF and PACF
#acf_plot <- ggAcf(residuals, lag.max = 30) + ggtitle("ACF of residuals")
#pacf_plot <- ggPacf(residuals, lag.max = 30) + ggtitle("PACF of residuals")

# QQ-plot
qq_plot <- ggplot(data.frame(resid = residuals), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "blue") + theme_minimal() +
  labs(title = "QQ-plot of residuals")

# Show all plots
grid.arrange(p_res, qq_plot, nrow = 2)
#grid.arrange(acf_plot, pacf_plot, nrow = 2)

#Compute AIC and BIC
# Log-likelihood value (positive because our function returns it this way)
logLik_val <- kf_logLik_dt(result$par, df)

# Number of parameters estimated = 7
k <- length(result$par)
n <- nrow(df)

# AIC and BIC
AIC_val <- -2 * logLik_val + 2 * k
BIC_val <- -2 * logLik_val + log(n) * k

cat("AIC:", AIC_val, "\n")
cat("BIC:", BIC_val, "\n")


