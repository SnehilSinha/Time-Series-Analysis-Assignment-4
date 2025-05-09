# Load libraries
library(tidyverse)
library(lubridate)
library(patchwork)
library(forecast)
library(gridExtra)

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

# Prepare data for Kalman filter
df <- data %>%
  transmute(
    Y  = Y_t,
    Ta = T_a_t,
    S  = Phi_s_t,
    I  = Phi_l_t
  )

# --- Kalman filter log-likelihood ---
kf_logLik_dt <- function(par, df) {
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  Sigma1lt <- matrix(exp(par[5]), 1, 1)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  C   <- matrix(1, 1, 1)
  Sigma2 <- matrix(exp(par[6]), 1, 1)
  X0  <- matrix(par[7], 1, 1)

  Y <- as.matrix(df[, "Y"])
  U <- as.matrix(df[, c("Ta", "S", "I")])
  Tn <- nrow(df)
  n <- nrow(A)

  x_est <- X0
  P_est <- diag(1e1, n)
  logLik <- 0

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t, ] - y_pred

    logLik <- logLik - 0.5 * (log(2 * pi) + log(det(S_t)) + t(innov) %*% solve(S_t, innov))

    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
  }

  as.numeric(logLik)
}

# --- Estimation function ---
estimate_dt <- function(start_par, df, lower = NULL, upper = NULL) {
  negLL <- function(par) { -kf_logLik_dt(par, df) }
  optim(
    par = start_par,
    fn = negLL,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 1000, trace = 1)
  )
}

# --- Kalman filter to return residuals and predictions ---
run_kalman_filter <- function(par, df) {
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  Sigma1lt <- matrix(exp(par[5]), 1, 1)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  C   <- matrix(1, 1, 1)
  Sigma2 <- matrix(exp(par[6]), 1, 1)
  X0  <- matrix(par[7], 1, 1)

  Y <- as.matrix(df[, "Y"])
  U <- as.matrix(df[, c("Ta", "S", "I")])
  Tn <- nrow(df)

  x_est <- X0
  P_est <- diag(1e1, 1)
  predictions <- numeric(Tn)
  residuals   <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t, ] - y_pred

    predictions[t] <- y_pred
    residuals[t] <- innov

    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(1) - K_t %*% C) %*% P_pred
  }

  list(predictions = predictions, residuals = residuals)
}

# --- Estimate parameters ---
start_par <- c(0.5, 0.1, 0.1, 0.1, log(0.5), log(0.5), 20)
lower     <- c(-1, -1, -1, -1, log(1e-6), log(1e-6), -100)
upper     <- c(1, 1, 1, 1, log(100), log(100), 100)

result <- estimate_dt(start_par, df, lower, upper)
est_par <- result$par
print(est_par)

# --- Run Kalman filter again to extract residuals and predictions ---
kf_result <- run_kalman_filter(est_par, df)
residuals <- kf_result$residuals
predicted <- kf_result$predictions

# --- Compute AIC and BIC ---
logLik_val <- kf_logLik_dt(est_par, df)
k <- length(est_par)
n <- nrow(df)

AIC_val <- -2 * logLik_val + 2 * k
BIC_val <- -2 * logLik_val + log(n) * k

cat("AIC:", AIC_val, "\n")
cat("BIC:", BIC_val, "\n")

# --- Plot residual diagnostics ---
# Residual time series
p_res <- ggplot(data.frame(t = 1:n, resid = residuals), aes(x = t, y = resid)) +
  geom_line(color = "darkred") +
  labs(title = "Residuals", y = "Residual", x = "Time")

# ACF and PACF
acf_plot <- ggAcf(residuals, lag.max = 30) + ggtitle("ACF of residuals")
pacf_plot <- ggPacf(residuals, lag.max = 30) + ggtitle("PACF of residuals")

# QQ-plot
qq_plot <- ggplot(data.frame(resid = residuals), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "blue") +
  labs(title = "QQ-plot of residuals")

# Observed vs predicted
obs_pred_plot <- ggplot(data.frame(
  time = 1:n,
  observed = df$Y,
  predicted = predicted
)) +
  geom_line(aes(x = time, y = observed), color = "black") +
  geom_line(aes(x = time, y = predicted), color = "blue") +
  labs(title = "Observed vs Predicted Temperature", y = "Y (°C)", x = "Time Step")

# Show all plots
print(p_res)
grid.arrange(p_res, qq_plot, acf_plot, pacf_plot, ncol = 2)
