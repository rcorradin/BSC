# This code was initially developed by Tommaso Rigon,
# the original version is available at 
# https://tommasorigon.github.io/CompStat/

library(coda)
library(survival)
library(tictoc) 

# AR example -------------------------------------------------------------------
# ------------------------------------------------------------------------------

R <- 300

par(mfrow = c(2,2))
# Stationary process
set.seed(123)
rho <- 0
y_0 <- numeric(R + 1)
y_0[1] <- rnorm(1, 30, 1)
for(r in 1:R){
  y_0[r + 1] = rho * y_0[r] + rnorm(1)
}
plot(y_0, type = "l")

# Stationary process
set.seed(123)
rho <- 0.5
y_05 <- numeric(R + 1)
y_05[1] <- rnorm(1, 30, 1)
for(r in 1:R){
  y_05[r + 1] = rho * y_05[r] + rnorm(1)
}
plot(y_05, type = "l")

# Non-stationary process
set.seed(123)
rho <- 1
y_1 <- numeric(R + 1)
y_1[1] <- rnorm(1, 30, 1)
for(r in 1:R){
  y_1[r + 1] = rho * y_1[r] + rnorm(1)
}
plot(y_1, type = "l")

# Non-stationary process
set.seed(123)
rho <- 1.01
y_101 <- numeric(R + 1)
y_101[1] <- rnorm(1, 30, 1)
for(r in 1:R){
  y_101[r + 1] = rho * y_101[r] + rnorm(1)
}
plot(y_101, type = "l")


# MH toy Gaussian --------------------------------------------------------------
# ------------------------------------------------------------------------------

norm_mcmc <- function(R, mu, sig, ep, x0) {
  # Initialization
  out <- numeric(R + 1)
  out[1] <- x0
  # Beginning of the chain
  x <- x0
  # Metropolis algorithm
  for(r in 1:R){
    # Proposed values
    xs <- x + runif(1,
                    -ep, ep)
    # Acceptance probability
    alpha <- min(dnorm(xs, mu, sig) / dnorm(x, mu, sig), 1)
    # Acceptance / rejection step
    accept <- rbinom(1, size = 1, prob = alpha)
    if(accept == 1) {
      x <- xs
    }
    out[r + 1] <- x
  }
  out
}

set.seed(123)
sim1 <- norm_mcmc(1000, mu = 2, sig = 5, ep = 100, x0 = 50)
sim2 <- norm_mcmc(1000, mu = 2, sig = 5, ep = 50, x0 = 50)
sim3 <- norm_mcmc(1000, mu = 2, sig = 5, ep = 10, x0 = 50)
sim4 <- norm_mcmc(1000, mu = 2, sig = 5, ep = 1, x0 = 50)

par(mfrow = c(2,2))
plot(sim1, type = "l", main = "ep = 100", ylab = "y", xlab = "iteration")
plot(sim2, type = "l", main = "ep = 50", ylab = "y", xlab = "iteration")
plot(sim3, type = "l", main = "ep = 10", ylab = "y", xlab = "iteration")
plot(sim4, type = "l", main = "ep = 1", ylab = "y", xlab = "iteration")

par(mfrow = c(2, 2))
acf(sim1)
acf(sim2)
acf(sim3)
acf(sim4)

# Simulate the MH chain
sim <- norm_mcmc(50000, mu = 2, sig = 5, ep = 10, x0 = 50)

# Identify a burn-in period
burn_in <- 1:200
sim <- sim[-c(burn_in)]
# Plot the results
par(mfrow = c(1,1))
hist(sim, breaks = 100, freq = FALSE)
curve(dnorm(x, 2, 5), add = T) # This is usually not known!

# HYBRID MH toy Gaussian -------------------------------------------------------
# ------------------------------------------------------------------------------

dbvnorm <- function(x, rho) {
  exp(-(x[1]^2 - 2 * rho * x[1] * x[2] + x[2]^2) / (2 * (1 - rho^2)))
}

bvnorm_mcmc <- function(R, rho, ep, x0) {
  out <- matrix(0, R + 1, 2)
  out[1, ] <- x0
  x <- x0
  for(r in 1:R){
    for(j in 1:2){
      xs <- x
      xs[j] <- x[j] + runif(1, -ep[j], ep[j])
      alpha <- min(dbvnorm(xs, rho) / dbvnorm(x, rho), 1) # Acceptance probability
      accept <- rbinom(1, size = 1, prob = alpha) # Acceptance / rejection step
      if(accept == 1) {
        x[j] <- xs[j]
      }
    }
    out[r + 1, ] <- x
  }
  out 
}

set.seed(123)
sim <- bvnorm_mcmc(1000, rho = 0.85, ep = c(2, 2), x0 = c(10, 10))
plot(sim, type = "b")

# GIBBS SAMPLER toy Gaussian ---------------------------------------------------
# ------------------------------------------------------------------------------

gibbs_R <- function(x, m_0, lambda_02, a_0, b_0, R, burn_in) {
  # Initialization
  n <- length(x); xbar <- mean(x)
  out <- matrix(0, R, 2)
  # Initial values for mu and sigma
  sigma2 <- var(x); mu <- xbar
  for (r in 1:(burn_in + R)) {
    # Sample mu
    sigma2_n <- 1 / (1 / lambda_02 + n / sigma2)
    mu_n <- sigma2_n * (m_0 / lambda_02 + n / sigma2 * xbar)
    mu <- rnorm(1, mu_n, sqrt(sigma2_n))
    # Sample sigma2
    a_n <- a_0 + 0.5 * n
    b_n <- b_0 + 0.5 * sum((x - mu)^2)
    sigma2 <- 1 / rgamma(1, a_n, b_n)
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- c(mu, sigma2)
    }
  }
  out
}

set.seed(123)
x <- rnorm(n = 25, mean = 1, sd = 1)
sim <- gibbs_R(x, m_0 = 0, lambda_02 = 100, a_0 = 3, b_0 = 10, R = 10000, burn_in = 2000)

fit_MCMC <- as.mcmc(sim) # Convert the matrix into a "coda" object
plot(fit_MCMC)

# SURVIVAL ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

t <- stanford2$time # Survival times
d <- stanford2$status # Censorship indicator

loglik <- function(t, d, alpha, beta) {
  log_hazard <- sum(d * ((alpha - 1) * log(t / beta) + log(alpha / beta)))
  log_survival <- sum(-(t / beta)^alpha)
  log_hazard + log_survival
}

logprior <- function(theta) {
  sum(dnorm(theta, 0, sqrt(100), log = TRUE))
}

logpost <- function(t, d, theta) {
  loglik(t, d, exp(theta[1]), exp(theta[2])) + logprior(theta)
}

# R represents the number of samples
# burn_in is the number of discarded samples
RMH <- function(R, burn_in, t, d) {
  out <- matrix(0, R, 2) # Initialize an empty matrix to store the values
  theta <- c(0, 0) # Initial values
  logp <- logpost(t, d, theta) # Log-posterior
  
  for (r in 1:(burn_in + R)) {
    theta_new <- rnorm(2, mean = theta, sd = 0.25) # Propose a new value
    logp_new <- logpost(t, d, theta_new)
    alpha <- min(1, exp(logp_new - logp))
    if (runif(1) < alpha) {
      theta <- theta_new # Accept the value
      logp <- logp_new
    }
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- theta
    }
  }
  out
}

R <- 50000
burn_in <- 5000

set.seed(123)
tic()
fit_MCMC <- RMH(R, burn_in, t, d)
toc()

library(coda)
fit_MCMC <- exp(fit_MCMC) # Back to the original parametrization
fit_MCMC <- as.mcmc(fit_MCMC) # Convert the matrix into a "coda" object
plot(fit_MCMC)
effectiveSize(fit_MCMC) # Effective sample size
1 - rejectionRate(fit_MCMC) # Acceptance rate

# Grid of values on which the survival function is computed
grid <- seq(0, 3700, length = 50)

# Initialized the empty vectors
S_mean <- numeric(length(grid))
S_upper <- numeric(length(grid))
S_lower <- numeric(length(grid))

for (i in 1:length(grid)) {
  S_mean[i] <- mean(pweibull(grid[i], shape = fit_MCMC[, 1], fit_MCMC[, 2], lower.tail = FALSE))
  S_lower[i] <- quantile(pweibull(grid[i], shape = fit_MCMC[, 1], fit_MCMC[, 2], lower.tail = FALSE), 0.025)
  S_upper[i] <- quantile(pweibull(grid[i], shape = fit_MCMC[, 1], fit_MCMC[, 2], lower.tail = FALSE), 0.975)
}

library(ggplot2)
data_plot <- data.frame(Time = grid, Mean = S_mean, Upper = S_upper, Lower = S_lower)
ggplot(data = data_plot, aes(x = Time, y = Mean, ymin = Lower, ymax = Upper)) +
  geom_line() +
  theme_bw() +
  ylab("Survival function") +
  ylim(c(0, 1)) +
  geom_ribbon(alpha = 0.1)
