library(ggplot2)
library(mvtnorm)

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


# IMPORTANCE SAMPLER -----------------------------------------------------------
# ------------------------------------------------------------------------------

R <- 30000 # Number of retained samples

# R represents the number of samples
# S is the covariance matrix of the multivariate Gaussian proposal

IS <- function(R, y, X, S, mu) {
  p <- ncol(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  proposed_values <- rmvnorm(R, mean = mu, sigma = S)
  weights <- apply(proposed_values, 1, function(x) logpost(x, y, X)) - 
    - dmvnorm(proposed_values, , mean = summary(fit_logit)$coefficients[,1], sigma = S, log = TRUE)
  list(values = proposed_values, weights = weights - logsumexp(weights))
}

# large proposal
S <- diag(100, ncol(X))
mu <- rep(0, 8)

set.seed(123)
IS_not_informed <- IS(R, y, X, S, mu)
est_coef <- colSums(IS_not_informed$values * exp(IS_small$weights))
idx <- sample(1:R, R, TRUE, exp(IS_not_informed$weights))
length(unique(idx))
ESS <- 1 / sum(exp(IS_not_informed$weights)^2)

# informed
S <- vcov(fit_logit)
mu <- coefficients(fit_logit)

set.seed(123)
IS_informed <- IS(R, y, X, S, mu)
est_coef <- colSums(IS_informed$values * exp(IS_small$weights))
idx <- sample(1:R, R, TRUE, exp(IS_informed$weights))
length(unique(idx))
ESS <- 1 / sum(exp(IS_informed$weights)^2)

# IMPORTANCE RESAMPLING --------------------------------------------------------
# ------------------------------------------------------------------------------

R <- 20
rho <- 0.5
tau <- 0.3
set.seed(123)

# Strategy 1
IR1 <- function(R, rho, tau){
  sample0 <- rnorm(R)
  weight0 <- dnorm(sample0) * ifelse(abs(sample0) <= tau, 1, 0)
  sample1 <- rho * sample0 + sqrt(1 - rho^2) * rnorm(R)
  list(sample0, sample1, weight0)
}

sampleIR1 <- IR1(R, rho, tau)

p1 <- ggplot(data.frame(x = sampleIR1[[1]], y = sampleIR1[[2]], q = rep(0, R), 
                  z = rep(1, R), fill = as.factor(ifelse(abs(sampleIR1[[1]]) <= tau, 1, 0)))) +
  geom_segment(aes(x = q, y = x, xend = z, yend = y), colour = "gray") +
  geom_point(aes(x = q, y = x, colour = fill), size = 4) +
  geom_point(aes(x = q, y = x), shape = 1, size = 4, colour = "black") +
  geom_point(aes(x = z, y = y, colour = fill), size = 4) +
  geom_point(aes(x = z, y = y), shape = 1, size = 4, colour = "black") +
  theme_bw() +
  theme(legend.position = "null") +
  xlab("time") +
  ylab("value")

# Strategy 1
IR2 <- function(R, rho, tau){
  sample0 <- rnorm(R)
  weight0 <- dnorm(sample0) * ifelse(abs(sample0) <= tau, 1, 0)
  temp_idx <- sample(1:R, R, T, weight0)
  sample1 <- rho * sample0[temp_idx] + sqrt(1 - rho^2) * rnorm(R)
  list(sample0, sample1, weight0, temp_idx)
}

sampleIR2 <- IR2(R, rho, tau)
p2 <- ggplot(data.frame(x = sampleIR2[[1]], y = sampleIR2[[2]], q = rep(0, R), from = sampleIR2[[1]][sampleIR2[[4]]], 
                  z = rep(1, R), fill = as.factor(ifelse(abs(sampleIR2[[1]]) <= tau, 1, 0)))) +
  geom_segment(aes(x = q, y = from, xend = z, yend = y), colour = "gray") +
  geom_point(aes(x = q, y = x, colour = fill), size = 4) +
  geom_point(aes(x = q, y = x), shape = 1, size = 4, colour = "black") +
  geom_point(aes(x = z, y = y), colour = "#00BFC4", size = 4) +
  geom_point(aes(x = z, y = y), shape = 1, size = 4, colour = "black") +
  theme_bw() +
  theme(legend.position = "null") +
  xlab("time") +
  ylab("value")

pdf(file = "plot.pdf", width = 7, height = 3.5)
ggpubr::ggarrange(p1, p2, ncol = 2)
dev.off()


# PIMA MALA RESAMPLER -----------------------------------------------------------
# ------------------------------------------------------------------------------

R <- 30000 # Number of retained samples

# R represents the number of samples
# S is the covariance matrix of the multivariate Gaussian proposal

IrS <- function(R, y, X, S, mu, sigma2) {
  p <- ncol(X)
  sigma <- sqrt(sigma2)
  A <- chol(S) 
  
  proposed_values0 <- rmvnorm(R, mean = mu, sigma = S)
  weights0 <- apply(proposed_values0, 1, function(x) logpost(x, y, X)) - 
    - dmvnorm(proposed_values0, mean = summary(fit_logit)$coefficients[,1], sigma = S, log = TRUE)
  
  proposed_values1 <- proposed_values0[sample(1:R, R, TRUE, exp(weights0 - logsumexp(weights0))),]
  proposed_values1 <- proposed_values1 + 
    t(apply(proposed_values1, 1, function(x) sigma2 / 2 * 
            c(S %*% lgradient(x, y, X)) + sigma * c(crossprod(A, rnorm(p)))))
  weights1 <- apply(proposed_values1, 1, function(x) logpost(x, y, X)) - 
    - dmvnorm(proposed_values1, mean = summary(fit_logit)$coefficients[,1], sigma = S, log = TRUE)
  
  list(values = proposed_values1, weights = weights1 - logsumexp(weights1))
}


# informed
S <- vcov(fit_logit)
mu <- coefficients(fit_logit)
sigma2 <- 0.1

set.seed(123)
IrS_not_informed <- IrS(R, y, X, S, mu, sigma2)
est_coef <- colSums(IrS_not_informed$values * exp(IrS_not_informed$weights))
idx <- sample(1:R, R, TRUE, exp(IrS_not_informed$weights))
length(unique(idx))
ESS <- 1 / sum(exp(IrS_not_informed$weights)^2)
