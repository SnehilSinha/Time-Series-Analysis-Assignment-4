# === Load libraries ===
library(tidyverse)
library(lubridate)
library(gridExtra)
library(ggplot2)
library(forecast)

# === Read and prepare the data ===
setwd("C:/Users/snehi/Documents/Time Series Analysis/Assignment4")
data <- read_csv("transformer_data.csv")

# Rename and prepare for modeling
data <- data %>%
  rename(Y_t = Y, T_a_t = Ta, Phi_s_t = S, Phi_l_t = I)

df <- data %>%
  transmute(
    Y  = Y_t,
    Ta = T_a_t,
    S  = Phi_s_t,
    I  = Phi_l_t
  )

# === 2D Kalman filter function ===
run_kalman_filter_2D <- function(par, df) {
  A <- matrix(par[1:4], nrow = 2)
  B <- matrix(par[5:10], nrow = 2)
  Sigma1_lt <- diag(exp(par[11:12]))
  Sigma1 <- Sigma1_lt %*% t(Sigma1_lt)
  C <- matrix(c(1, 0), nrow = 1)
  Sigma2 <- matrix(exp(par[13]), 1, 1)
  X0 <- matrix(par[14:15], nrow = 2)

  Y <- as.matrix(df$Y)
  U <- as.matrix(df[, c("Ta", "S", "I")])
  Tn <- nrow(df)

  n <- nrow(A)
  x_est <- X0
  P_est <- diag(10, n)
  predictions <- numeric(Tn)
  residuals <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred

    predictions[t] <- y_pred
    residuals[t] <- innov

    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
  }

  list(predictions = predictions, residuals = residuals)
}

# === Log-likelihood function ===
kf_logLik_dt <- function(par, df) {
  A <- matrix(par[1:4], nrow = 2)
  B <- matrix(par[5:10], nrow = 2)
  Sigma1_lt <- diag(exp(par[11:12]))
  Sigma1 <- Sigma1_lt %*% t(Sigma1_lt)
  C <- matrix(c(1, 0), nrow = 1)
  Sigma2 <- matrix(exp(par[13]), 1, 1)
  X0 <- matrix(par[14:15], nrow = 2)

  Y <- as.matrix(df$Y)
  U <- as.matrix(df[, c("Ta", "S", "I")])
  Tn <- nrow(df)

  n <- nrow(A)
  x_est <- X0
  P_est <- diag(10, n)
  logLik <- 0

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred

    logLik <- logLik - 0.5 * (log(det(S_t)) + t(innov) %*% solve(S_t) %*% innov)

    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
  }

  return(as.numeric(logLik))
}

# === Estimation wrapper ===
estimate_dt <- function(start_par, df, lower = NULL, upper = NULL) {
  negLL <- function(par) { -kf_logLik_dt(par, df) }
  optim(
    par     = start_par,
    fn      = negLL,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = 1000, trace = 1)
  )
}

# === Initial parameters (15) ===
start_par <- c(
  0.9, 0.1,  # A[1,]
  0.1, 0.9,  # A[2,]
  0.1, 0.1, 0.1,  # B[1,]
  0.1, 0.1, 0.1,  # B[2,]
  log(0.1), log(0.1),  # state noise std devs (log)
  log(0.5),            # obs noise std dev (log)
  20, 20               # X0
)

lower <- c(rep(-2, 10), rep(log(1e-6), 3), rep(-100, 2))
upper <- c(rep(2, 10),  rep(log(100),  3), rep(100, 2))

# === Estimate model parameters ===
result_2D <- estimate_dt(start_par, df, lower, upper)

# === Retry if optimizer didn’t converge ===
if (result_2D$convergence != 0) {
  cat("Re-running optimization due to non-convergence...\n")
  result_2D <- estimate_dt(result_2D$par, df, lower, upper)
}

print(result_2D$par)

# === Run Kalman filter with estimated parameters ===
kf_output <- run_kalman_filter_2D(result_2D$par, df)
residuals <- kf_output$residuals

# === Residual plots ===
# Residual time series
p_res <- ggplot(data.frame(t = 1:length(residuals), resid = residuals), aes(x = t, y = resid)) +
  geom_line(color = "darkred") + theme_minimal() +
  labs(title = "Residuals", y = "Residual", x = "Time")

# QQ-plot
qq_plot <- ggplot(data.frame(resid = residuals), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "blue") + theme_minimal() +
  labs(title = "QQ-plot of Residuals")


# ACF & PACF plots
acf_plot <- ggAcf(residuals, lag.max = 30) + ggtitle("ACF of Residuals")
pacf_plot <- ggPacf(residuals, lag.max = 30) + ggtitle("PACF of Residuals")

#plot all 4
grid.arrange(p_res, qq_plot, acf_plot, pacf_plot, ncol = 2)


# === AIC and BIC ===
logLik_val <- kf_logLik_dt(result_2D$par, df)
k <- length(result_2D$par)
n <- nrow(df)
AIC_val <- -2 * logLik_val + 2 * k
BIC_val <- -2 * logLik_val + log(n) * k

cat("AIC:", AIC_val, "\n")
cat("BIC:", BIC_val, "\n")
