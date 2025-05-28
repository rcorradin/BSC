library(ggplot2)
library(mvtnorm)
library(coda)

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

# Logit model
fit_logit <- glm(type ~ X - 1, family = binomial(link = "logit"), data = Pima)
summary(fit_logit)

# LOGSUMEXP --------------------------------------------------------------------
# ------------------------------------------------------------------------------

logsumexp <- function(v){
  M <- max(v)
  out <- M + log(sum(exp(v - M)))
  out
}

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


# CAVI ALGORITHM ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# Compute the log-determinant of a matrix
ldet <- function(X) {
  if (!is.matrix(X)) {
    return(log(X))
  }
  determinant(X, logarithm = TRUE)$modulus
}

logit_CAVI <- function(y, X, B, b, tol = 1e-16, maxiter = 10000) {
  lowerbound <- numeric(maxiter)
  p <- ncol(X)
  n <- nrow(X)
  
  P <- solve(B)
  Pb <- c(P %*% b)
  Pdet <- ldet(P)
  
  # Initialization for omega equal to 0.25
  P_vb <- crossprod(X * rep(1 / 4, n), X) + P
  Sigma_vb <- solve(P_vb)
  mu_vb <- Sigma_vb %*% (crossprod(X, y - 0.5) + Pb)
  eta <- c(X %*% mu_vb)
  xi <- sqrt(eta^2 + rowSums(X %*% Sigma_vb * X))
  omega <- tanh(xi / 2) / (2 * xi)
  omega[is.nan(omega)] <- 0.25
  
  lowerbound[1] <- 0.5 * p + 0.5 * ldet(Sigma_vb) + 0.5 * Pdet - 0.5 * t(mu_vb - b) %*% P %*% (mu_vb - b) + sum((y - 0.5) * eta + log(plogis(xi)) - 0.5 * xi) - 0.5 * sum(diag(P %*% Sigma_vb))
  
  # Iterative procedure
  for (t in 2:maxiter) {
    P_vb <- crossprod(X * omega, X) + P
    Sigma_vb <- solve(P_vb)
    mu_vb <- Sigma_vb %*% (crossprod(X, y - 0.5) + Pb)
    
    # Update of xi
    eta <- c(X %*% mu_vb)
    xi <- sqrt(eta^2 + rowSums(X %*% Sigma_vb * X))
    omega <- tanh(xi / 2) / (2 * xi)
    omega[is.nan(omega)] <- 0.25
    
    lowerbound[t] <- 0.5 * p + 0.5 * ldet(Sigma_vb) + 0.5 * Pdet - 0.5 * t(mu_vb - b) %*% P %*% (mu_vb - b) + sum((y - 0.5) * eta + log(plogis(xi)) - 0.5 * xi) - 0.5 * sum(diag(P %*% Sigma_vb))
    
    if (abs(lowerbound[t] - lowerbound[t - 1]) < tol) {
      return(list(
        mu = c(mu_vb), Sigma = matrix(Sigma_vb, p, p),
        Convergence = cbind(
          Iteration = (1:t) - 1,
          Lowerbound = lowerbound[1:t]
        ), xi = xi
      ))
    }
  }
  stop("The algorithm has not reached convergence")
}

b <- rep(0, 8)
B <- diag(100, 8)
fit_CAVI <- logit_CAVI(y, X, B, b)
est_coef <- fit_CAVI$mu

# ABC REJECTION ----------------------------------------------------------------
# ------------------------------------------------------------------------------

ABCrejection <- function(R, lambda, y, X, S, mu){
  p <- ncol(X)
  out <- matrix(0, R, p) 
  A <- chol(S) 
  C <- optim(f = function(x) (dmvnorm(x, mean = rep(0, 8), sigma = diag(100, 8), log = T) - 
               dmvnorm(x, mean = mu, sigma = S, log = T)), par = rep(0, 8), method = "L-BFGS-B",
               lower = rep(-100, 8), upper = rep(100, 8))$value
  trials <- 0
  
  for(r in 1:R){
    samp <- TRUE
    while(samp){
      trials <- trials + 1
      temp_param <- as.vector(rnorm(8) %*% A + mu)
      temp_pred <- X %*% temp_param
      temp_probs <- exp(temp_pred) / (1 + exp(temp_pred))
      s <- sapply(temp_probs, function(x) rbinom(1, c(1, 0), c(x, 1 - x)))
      weight <- mean(y == s)^lambda * 
        exp(dmvnorm(temp_param, mean = rep(0, 8), sigma = diag(100, 8), log = T) -
        dmvnorm(temp_param, mean = mu, sigma = S, log = T) - C)
      if(runif(1) < weight){
        out[r,] <- temp_param
        samp <- FALSE
      }
    }
  }
  list(param = out, acc = R / trials)
}

# not informed

R <- 3000
S <- diag(100, 8)
mu <- rep(0, 8)
lambda <- 3

set.seed(123)
ABCrejectionNI <- ABCrejection(R, lambda, y, X, S, mu)
colMeans(ABCrejectionNI$param)
ABCrejectionNI$acc

# informed

R <- 3000
S <- vcov(fit_logit)
mu <- coefficients(fit_logit)
lambda <- 3

set.seed(123)
ABCrejectionI <- ABCrejection(R, lambda, y, X, S, mu)
colMeans(ABCrejectionI$param)
ABCrejectionI$acc


# ABC IMPORTANCE ---------------------------------------------------------------
# ------------------------------------------------------------------------------

ABCimportance <- function(R, lambda, y, X, S, mu){
  p <- ncol(X)
  out <- matrix(0, R, p)
  weights <- rep(0, R)
  A <- chol(S) 
  
  for(r in 1:R){
    out[r,] <- t(rnorm(8) %*% A + mu)
    temp_pred <- X %*% out[r,]
    temp_probs <- exp(temp_pred) / (1 + exp(temp_pred))
    s <- sapply(temp_probs, function(x) rbinom(1, c(1, 0), c(x, 1 - x)))
    weights[r] <- mean(y == s)^lambda * 
      dmvnorm(out[r,], rep(0,8), diag(100, 8)) / dmvnorm(out[r,], mu, S)
  }
  list(param = out, weights = weights / sum(weights))
}

# not informed

R <- 3000
S <- diag(100, 8)
mu <- rep(0, 8)
lambda <- 3

set.seed(123)
ABCimportanceNI <- ABCimportance(R, lambda, y, X, S, mu)
est_param <- colSums(ABCimportanceNI$param * ABCimportanceNI$weights)
ESS <-  1 / sum(ABCimportanceNI$weights^2)

# informed

R <- 3000
S <- vcov(fit_logit)
mu <- coefficients(fit_logit)
lambda <- 3

set.seed(123)
ABCimportanceI <- ABCimportance(R, lambda, y, X, S, mu)
est_param <- colSums(ABCimportanceI$param * ABCimportanceI$weights)
ESS <-  1 / sum(ABCimportanceI$weights^2)
