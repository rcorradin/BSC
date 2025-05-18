# This code was initially developed by Tommaso Rigon,
# the original version is available at 
# https://tommasorigon.github.io/CompStat/

library(coda)
library(MASS)
library(Rcpp)

Pima <- rbind(MASS::Pima.tr, MASS::Pima.te)
head(Pima)

# A population of women who were at least 21 years old, of Pima Indian heritage, 
# and living near Phoenix, Arizona, was tested for diabetes according to World 
# Health Organization criteria. The data were collected by the US National 
# Institute of Diabetes and Digestive and Kidney Diseases. We used the 532 
# complete records after dropping the (mainly missing) data on serum insulin.

y <- as.numeric(Pima$type == "Yes") # Binary outcome
X <- model.matrix(type ~ . - 1, data = Pima) # Design matrix
X <- cbind(1, scale(X)) # Standardize the design matrix, add the intercept
head(X)

# Probit model
fit_probit <- glm(type ~ X - 1, family = binomial(link = "probit"), data = Pima)
summary(fit_probit)
# Logit model
fit_logit <- glm(type ~ X - 1, family = binomial(link = "logit"), data = Pima)
summary(fit_logit)

# GOING BAYESIAN ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Some quantities needed in the sample

# Loglikelihood of a logistic regression model
loglik <- function(beta, y, X) {
  eta <- c(X %*% beta)
  sum(y * eta - log(1 + exp(eta)))
}

# Logposterior
logpost <- function(beta, y, X) {
  loglik(beta, y, X) + sum(dnorm(beta, 0, 10, log = T))
}

# Gradient of the logposterior
lgradient <- function(beta, y, X) {
  probs <- plogis(c(X %*% beta))
  loglik_gr <- c(crossprod(X, y - probs))
  prior_gr <- -beta / 100
  loglik_gr + prior_gr
}

# OUTPUT -----------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Summary table
summary_tab <- matrix(0, nrow = 5, ncol = 4)
colnames(summary_tab) <- c("Seconds", "Average ESS", "Average ESS per sec", "Average acceptance rate")
rownames(summary_tab) <- c("Vanilla MH", "Laplace MH", "Adaptive MH", "Metropolis within Gibbs", "Adaptive Metropolis within Gibbs")

# VANILLA MH -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# R represents the number of samples
# burn_in is the number of discarded samples
# S is the covariance matrix of the multivariate Gaussian proposal
RMH <- function(R, burn_in, y, X, S) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- rep(0, p) # Initial values
  logp <- logpost(beta, y, X)
  
  # Eigen-decomposition
  eig <- eigen(S, symmetric = TRUE)
  A1 <- t(eig$vectors) * sqrt(eig$values)
  
  # Starting the Gibbs sampling
  for (r in 1:(burn_in + R)) {
    beta_new <- beta + c(matrix(rnorm(p), 1, p) %*% A1)
    logp_new <- logpost(beta_new, y, X)
    alpha <- min(1, exp(logp_new - logp))
    if (runif(1) < alpha) {
      logp <- logp_new
      beta <- beta_new # Accept the value
    }
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}

R <- 30000 # Number of retained samples
burn_in <- 30000 # Burn-in period
set.seed(123)

# Covariance matrix of the proposal
S <- diag(1e-3, ncol(X))

# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(RMH(R, burn_in, y, X, S)) # Convert the matrix into a "coda" object
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[1, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of a couple of parameters
plot(fit_MCMC[, 3:4])

# LAPLACE MH -------------------------------------------------------------------
# ------------------------------------------------------------------------------

set.seed(123)

# Running the MCMC
start.time <- Sys.time()

# Covariance matrix is selected via Laplace approximation
fit_logit <- glm(type ~ X - 1, family = binomial(link = "logit"), data = Pima)
p <- ncol(X)
S <- 2.38^2 * vcov(fit_logit) / p
# MCMC
fit_MCMC <- as.mcmc(RMH(R, burn_in, y, X, S)) # Convert the matrix into a "coda" object
end.time <- Sys.time()

time_in_sec <- as.numeric(end.time - start.time)

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[2, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# ADAPTIVE METROPOLIS ----------------------------------------------------------
# ------------------------------------------------------------------------------

# R represents the number of samples
# burn_in is the number of discarded samples
# S is the covariance matrix of the multivariate Gaussian proposal
RMH_Adaptive <- function(R, burn_in, y, X) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- rep(0, p) # Initial values
  logp <- logpost(beta, y, X)
  epsilon <- 1e-6 # Inital value for the covariance matrix
  
  # Initial matrix S
  S <- diag(epsilon, p)
  Sigma_r <- diag(0, p)
  mu_r <- beta
  
  for (r in 1:(burn_in + R)) {
    
    # Updating the covariance matrix
    if(r > 1){
      Sigma_r <- (r - 2) / (r - 1) * Sigma_r + tcrossprod(beta - mu_r) / r
      mu_r <- (r - 1) / r * mu_r + beta / r
      S <- 2.38^2 * Sigma_r / p + diag(epsilon, p)
    }
    
    # Eigen-decomposition
    eig <- eigen(S, symmetric = TRUE)
    A1 <- t(eig$vectors) * sqrt(eig$values)
    
    beta_new <- beta + c(matrix(rnorm(p), 1, p) %*% A1)
    logp_new <- logpost(beta_new, y, X)
    alpha <- min(1, exp(logp_new - logp))
    if (runif(1) < alpha) {
      logp <- logp_new
      beta <- beta_new # Accept the value
    }
    
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}

set.seed(123)

# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(RMH_Adaptive(R = R, burn_in = burn_in, y, X))
end.time <- Sys.time()

time_in_sec <- as.numeric(end.time - start.time)

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[3, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# METROPOLIS WITHIN GIBBS ------------------------------------------------------
# ------------------------------------------------------------------------------

RMH_Gibbs <- function(R, burn_in, y, X, se) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- beta_new <- rep(0, p) # Initial values
  logp <- logpost(beta, y, X)
  
  for (r in 1:(burn_in + R)) {
    for (j in 1:p) {
      beta_new[j] <- beta[j] + rnorm(1, 0, se[j])
      logp_new <- logpost(beta_new, y, X)
      alpha <- min(1, exp(logp_new - logp))
      if (runif(1) < alpha) {
        logp <- logp_new
        beta[j] <- beta_new[j] # Accept the value
      }
    }
    
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}

p <- ncol(X)
se <- sqrt(rep(1e-4, p))

set.seed(123)
# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(RMH_Gibbs(R = R, burn_in = burn_in, y, X, se)) # Convert the matrix into a "coda" object
end.time <- Sys.time()

time_in_sec <- as.numeric(end.time - start.time)

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[4, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# ADAPTIVE METROPOLIS WITHIN GIBBS ---------------------------------------------
# ------------------------------------------------------------------------------

RMH_Gibbs_Adaptive <- function(R, burn_in, y, X, target = 0.44) {
  
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- beta_new <- rep(0, p) # Initial values
  logp <- logpost(beta, y, X)
  
  epsilon <- rep(0, p) # Initial value
  accepted <- numeric(p) # Vector of accepted values
  batch <- 1
  
  for (r in 1:(burn_in + R)) {
    
    # Do we need to update the parameters?
    if (batch == 50) {
      for (j in 1:p) {
        # Adapting the standard errors
        if ((accepted[j] / 50) > target) {
          epsilon[j] <- epsilon[j] + min(0.01, sqrt(1 / r))
        } else {
          epsilon[j] <- epsilon[j] - min(0.01, sqrt(1 / r))
        }
      }
      # Restart the cycle - Erase everything
      accepted <- numeric(p) # Vector of accepted values
      batch <- 0
    }
    
    # Increment the batch
    batch <- batch + 1
    
    for (j in 1:p) {
      beta_new[j] <- beta[j] + rnorm(1, 0, exp(epsilon[j]))
      logp_new <- logpost(beta_new, y, X)
      alpha <- min(1, exp(logp_new - logp))
      if (runif(1) < alpha) {
        logp <- logp_new
        beta[j] <- beta_new[j] # Accept the value
        accepted[j] <- accepted[j] + 1
      }
    }
    
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}

set.seed(123)
# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(RMH_Gibbs_Adaptive(R = R, burn_in = burn_in, y, X)) # Convert the matrix into a "coda" object
end.time <- Sys.time()

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[5, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept

plot(fit_MCMC[, 1:2])

# ------------------------------------------------------------------------------
# Print all the summaries

print(summary_tab)

# DYNAMIC-BASED METHODS --------------------------------------------------------
# ------------------------------------------------------------------------------

# Summary table for the six considered methods
summary_tab <- matrix(0, nrow = 4, ncol = 4)
colnames(summary_tab) <- c("Seconds", "Average ESS", "Average ESS per sec", "Average acceptance rate")
rownames(summary_tab) <- c("MH Laplace + Rcpp", "MALA", "MALA tuned", "HMC")

# Metropolis-Hastings (Laplace) and Rcpp
# ------------------------------------------------------------------------------

sourceCpp("RWMH.cpp")
R <- 30000
burn_in <- 5000
set.seed(123)

# Covariance matrix is selected via Laplace approximation
fit_logit <- glm(y ~ X - 1, family = binomial(link = "logit"))
p <- ncol(X)
S <- 2.38^2 * vcov(fit_logit) / p

# Running the MCMC
start.time <- Sys.time()
# MCMC
fit_MCMC <- as.mcmc(RMH_arma(R, burn_in, y, X, S)) # Convert the matrix into a "coda" object
end.time <- Sys.time()

time_in_sec <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[1, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# MALA Algorithm ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# R represents the number of samples
# burn_in is the number of discarded samples
# epsilon, S are tuning parameter
MALA <- function(R, burn_in, y, X, epsilon, S) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- rep(0, p) # Initial values
  A <- chol(S) # Cholesky of S
  S1 <- solve(S) # Inverse of S
  
  lgrad <- c(S %*% lgradient(beta, y, X)) # Compute the gradient
  logp <- logpost(beta, y, X)
  
  sigma2 <- epsilon^2 / p^(1 / 3)
  sigma <- sqrt(sigma2)
  
  # Starting the Gibbs sampling
  for (r in 1:(burn_in + R)) {
    beta_new <- beta + sigma2 / 2 * lgrad + sigma * c(crossprod(A, rnorm(p)))
    
    logpnew <- logpost(beta_new, y, X)
    lgrad_new <- c(S %*% lgradient(beta_new, y, X))
    
    diffold <- beta - beta_new - sigma2 / 2 * lgrad_new
    diffnew <- beta_new - beta - sigma2 / 2 * lgrad
    
    qold <- -diffold %*% S1 %*% diffold / (2 * sigma2)
    qnew <- -diffnew %*% S1 %*% diffnew / (2 * sigma2)
    
    alpha <- min(1, exp(logpnew - logp + qold - qnew))
    if (runif(1) < alpha) {
      logp <- logpnew
      lgrad <- lgrad_new
      beta <- beta_new # Accept the value
    }
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}


R <- 30000
burn_in <- 5000
set.seed(123)
epsilon <- 0.0017 # After some trial and error

# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(MALA(R = R, burn_in = burn_in, y, X, epsilon, S = diag(ncol(X)))) # Convert the matrix into a "coda" object
end.time <- Sys.time()
time_in_sec <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[2, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# MALA algorithm with pre-conditioning -----------------------------------------
# ------------------------------------------------------------------------------

R <- 30000
burn_in <- 5000
set.seed(123)
epsilon <- 1.68 # After some trial and error

# Covariance matrix is selected via Laplace approximation
fit_logit <- glm(y ~ X - 1, family = binomial(link = "logit"))
S <- vcov(fit_logit)

# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(MALA(R = R, burn_in = burn_in, y, X, epsilon, S)) # Convert the matrix into a "coda" object
end.time <- Sys.time()
time_in_sec <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[3, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# Hamiltonian Monte Carlo ------------------------------------------------------
# ------------------------------------------------------------------------------

HMC <- function(R, burn_in, y, X, epsilon, S, L = 10) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  beta <- rep(0, p) # Initial values
  logp <- logpost(beta, y, X) # Initial log-posterior
  S1 <- solve(S)
  A1 <- chol(S1)
  
  # Starting the Gibbs sampling
  for (r in 1:(burn_in + R)) {
    P <- c(crossprod(A1, rnorm(p))) # Auxiliary variables
    logK <- c(P %*% S %*% P / 2) # Kinetic energy at the beginning of the trajectory
    
    # Make a half step for momentum at the beginning
    beta_new <- beta
    Pnew <- P + epsilon * lgradient(beta_new, y, X) / 2
    
    # Alternate full steps for position and momentum
    for (l in 1:L) {
      # Make a full step for the position
      beta_new <- beta_new + epsilon * c(S %*% Pnew)
      # Make a full step for the momentum, except at the end of the trajectory
      if (l != L) Pnew <- Pnew + epsilon * lgradient(beta_new, y, X)
    }
    # Make a half step for momentum at the end.
    Pnew <- Pnew + epsilon * lgradient(beta_new, y, X) / 2
    
    # Negate momentum at the end of the trajectory to make the proposal symmetric
    Pnew <- - Pnew
    
    # Evaluate potential and kinetic energies at the end of the trajectory
    logpnew <- logpost(beta_new, y, X)
    logKnew <- Pnew %*% S %*% Pnew / 2 
    
    # Accept or reject the state at the end of the trajectory, returning either
    # the position at the end of the trajectory or the initial position
    if (runif(1) < exp(logpnew - logp + logK - logKnew)) {
      logp <- logpnew
      beta <- beta_new # Accept the value
    }
    
    # Store the values after the burn-in period
    if (r > burn_in) {
      out[r - burn_in, ] <- beta
    }
  }
  out
}


R <- 30000
burn_in <- 5000
set.seed(123)
epsilon <- 0.05 # After some trial and error
L <- 30

# Covariance matrix is selected via Laplace approximation
fit_logit <- glm(y ~ X - 1, family = binomial(link = "logit"))
S <- vcov(fit_logit)

# Running the MCMC
start.time <- Sys.time()
fit_MCMC <- as.mcmc(HMC(R = R, burn_in = burn_in, y, X, epsilon, S, L)) # Convert the matrix into a "coda" object
end.time <- Sys.time()
time_in_sec <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Diagnostic
summary(effectiveSize(fit_MCMC)) # Effective sample size
summary(R / effectiveSize(fit_MCMC)) # Integrated autocorrelation time
summary(1 - rejectionRate(fit_MCMC)) # Acceptance rate

# Summary statistics
summary_tab[4, ] <- c(
  time_in_sec, mean(effectiveSize(fit_MCMC)),
  mean(effectiveSize(fit_MCMC)) / time_in_sec,
  1 - mean(rejectionRate(fit_MCMC))
)

# Traceplot of the intercept
plot(fit_MCMC[, 1:2])

# ------------------------------------------------------------------------------
# Print all the summaries

print(round(summary_tab, digits = 3))
