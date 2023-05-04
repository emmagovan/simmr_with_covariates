# A simple version of the covariate model adapted from
# 'simplest_covariate_model.Rmd'
#With two covariates

# Load in package
library(simmr)
library(readxl)
library(tidyverse)

# Source in the generic functions
source("FF_VB_generic_functions_correct.R")


# Extract data ------------------------------------------------------------

path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] |> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]

# Put them into the right names
n <- nrow(consumer)
n_isotopes <- 2
K <- nrow(sources)
mu_s <- sources[, c(2, 3)] #+ disc[, c(2,4)]
sigma_s <- sources[, c(4, 5)]
mu_c <- TEFs[, c(2, 3)]
sigma_c <- TEFs[, c(4, 5)]
q <- conc[, c(2:3)]
x1 <- consumer$Sex
x2 <- consumer$Age
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

# Variables for the FFVB
S <- 100
mu_alpha_zero <- c(rep(0, K))
mu_beta1_zero <- c(rep(0, K))
mu_beta2_zero <- c(rep(0, K))
sigma_alpha_zero <- c(rep(1, K))
sigma_beta1_zero <- c(rep(1, K))
sigma_beta2_zero <- c(rep(1, K))
n_isotopes <- ncol(mu_c)
c_0 <- c(rep(1, n_isotopes))
d_0 <- c(rep(1, n_isotopes))
lambda <- c(
  rep(0, K), #alpha
  rep(1, K * (K + 1) / 2), 
  rep(0, K), #beta_1
  rep(1, K * (K + 1) / 2),
  rep(0, K), #beta_2
  rep(1, K * (K + 1) / 2),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

# function to extract lambdas --------------------------------------------
lambda_extract <- function(n_covariates, K, n_isotopes){
  mat_size = K * (K+1) /2
  mu_beta = matrix(data = NA, nrow = n_covariates, ncol = K)
  sigma_beta = matrix(data = NA, nrow = n_covariates, ncol = mat_size)
  
  for(i in 1:n_covariates){
    mu_beta[i,] = (i * mat_size + i * K +1):(i * mat_size + i * K + K)
    sigma_beta[i,] = (i * mat_size + (i+1) *K +1): (i * mat_size + (i+1) *K +mat_size)
  }
  
  c = (sigma_beta[n_covariates, mat_size] + 1):(sigma_beta[n_covariates, mat_size] + n_isotopes)
  d = (sigma_beta[n_covariates, mat_size] + n_isotopes + 1):(sigma_beta[n_covariates, mat_size] + 2 * n_isotopes)
  
  return(list(mu_alpha = 1:K,
              sigma_alpha = (K+1):(K + mat_size),
              mu_beta = mu_beta,
              sigma_beta = sigma_beta,
              c = c,
              d = d
  ))
}
#just use here for now - then build in to r
lambda_index <- lambda_extract(2, K, n_isotopes)


# sim_theta ---------------------------------------------------------------

sim_theta <- function(S = 100, lambda) {
  # lambda contains the parameters
  # mean_alpha, 1:K
  # chol(Sigma_alpha), (K + 1):(K + (K * (K + 1)) / 2)
  # mean_beta_1, (K + (K * (K + 1) / 2) + 1):(2 * K + (K * (K + 1) / 2))
  # chol(Sigma_beta_1), (2 * K + (K * (K + 1) / 2) + 1):(K * (K + 1) + 2 * K)
  # mean_beta_2, (K * (K + 1) + 2 * K + 1):(K * (K + 1) + 3 * K)
  # chol(Sigma_beta), (K * (K + 1) + 3 * K + 1):(K * (K + 1) + 3 * K + K * (K + 1)/2)
  # sigma_j_shape, (K * (K + 1) + 3 * K + K * (K + 1)/2 +1):(K * (K + 1) + 3 * K + K * (K + 1)/2 +n_isotopes)
  # sigma_j_scale, (K * (K + 1) + 3 * K + K * (K + 1)/2 +n_isotopes +1):(K * (K + 1) + 3 * K + K * (K + 1)/2 + 2* n_isotopes)
  
  mean_alpha <- lambda[lambda_index$mu_alpha]
  
  # K*(K-1) precision terms
  chol_prec_alpha <- matrix(0, nrow = K, ncol = K)
  chol_prec_alpha[upper.tri(chol_prec_alpha, diag = TRUE)] <-
    lambda[lambda_index$sigma_alpha]
  
  mean_beta_1 <- lambda[lambda_index$mu_beta[1,]]
  chol_prec_beta_1 <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta_1[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
    lambda[lambda_index$sigma_beta[1,]]
  
  mean_beta_2 <- lambda[lambda_index$mu_beta[2,]]
  chol_prec_beta_2 <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta_2[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
    lambda[lambda_index$sigma_beta[2,]]
  
  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    t(rMVNormC(S, mu = mean_alpha, U = chol_prec_alpha)),
    t(rMVNormC(S, mu = mean_beta_1, U = chol_prec_beta_1)), 
    t(rMVNormC(S, mu = mean_beta_2, U = chol_prec_beta_2)), 
    matrix(
      rgamma(S * n_isotopes,
             shape = lambda[lambda_index$c],
             rate = lambda[lambda_index$d]
      ),
      nrow = S,
      ncol = n_isotopes,
      byrow = TRUE
    )
  )
  
  return(theta)
}
theta <- sim_theta(S, lambda)
# Theta is alpha (K of these), beta_1 (K of these), beta_2 (K of these), and sigma (J of these)
# lambda is mu_alpha, chol(sigma_alpha), mu_beta_1, chol(sigma_beta_1), mu_beta_2,
# chol(sigma_beta_2), c, d

# h -----------------------------------------------------------------------

# Log of likelihood added to prior
h <- function(theta) {
  # Create alpha, beta and sigma
  alpha <- theta[1:K]
  beta_1 <- theta[(K + 1):(2 * K)]
  beta_2 <- theta[(2 * K +1):(3 * K)]
  sigma <- theta[(3 * K + 1):(3 * K + 2)]
  f <- matrix(NA, ncol = K, nrow = n)
  for (i in 1:n) {
    for (k in 1:K) {
      f[i, k] <- alpha[k] + beta_1[k] * x1[i] + beta_2[k] * x2[i]
    }
  }
  p <- matrix(NA, ncol = K, nrow = n)
  for (i in 1:n) {
    p[i, ] <- exp(f[i, 1:K]) / (sum((exp(f[i, 1:K]))))
  }
  
  ## Im not sure this bit needs to be looped over?
  hold <- 0
  for (i in 1:n) {
    for (j in 1:n_isotopes) {
      hold <- hold + sum(dnorm(y[i, j],
                               mean = sum(p[i, ] * q[, j] * (mu_s[, j] + mu_c[, j])) / 
                                 sum(p[i, ] * q[, j]),
                               sd = sqrt(sum(p[i, ]^2 * q[, j]^2 * (sigma_s[, j]^2 + sigma_c[, j]^2)) / 
                                           sum(p[i, ]^2 * q[, j]^2) + sigma[j]^2),
                               log = TRUE
      ))
    }
  }
  
  return(hold + sum(dnorm(alpha, mu_alpha_zero, sigma_alpha_zero, log = TRUE)) +
           sum(dnorm(beta_1, mu_beta1_zero, sigma_beta1_zero, log = TRUE)) +
           sum(dnorm(beta_2, mu_beta2_zero, sigma_beta2_zero, log = TRUE)) +
           sum(dgamma(sigma, shape = c_0, rate = d_0, log = TRUE)))
}
h(theta[1, ])


# log_q -------------------------------------------------------------------

log_q <- function(lambda, theta) {
  mean_alpha <- lambda[lambda_index$mu_alpha]
  chol_prec_alpha <- matrix(0, nrow = K, ncol = K)
  chol_prec_alpha[upper.tri(chol_prec_alpha, diag = TRUE)] <-
    lambda[lambda_index$sigma_alpha]
  
  mean_beta_1 <- lambda[lambda_index$mu_beta[1,]]
  chol_prec_beta_1 <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta_1[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
    lambda[lambda_index$sigma_beta[1,]]
  
  mean_beta_2 <- lambda[lambda_index$mu_beta[2,]]
  chol_prec_beta_2 <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta_2[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
    lambda[lambda_index$sigma_beta[2,]]
  
  shape_sigma <- lambda[lambda_index$c]
  rate_sigma <- lambda[lambda_index$d]
  
  # Extract alpha, beta and sigma from theta
  alpha <- theta[1:K]
  beta_1 <- theta[(K + 1):(2 * K)]
  beta_2 <- theta[(2 * K + 1):(3 * K)]
  sigma <- theta[(3 * K + 1):(3 * K + n_isotopes)]
  
  # This is a placeholder for more advanced code using chol_prec directly
  # prec <- crossprod(chol_prec)
  p1 <- matrix(alpha - mean_alpha, nrow = 1) %*% t(chol_prec_alpha)
  p2 <- matrix(beta_1 - mean_beta_1, nrow = 1) %*% t(chol_prec_beta_1)
  p3 <- matrix(beta_2 - mean_beta_2, nrow = 1) %*% t(chol_prec_beta_2)
  # log_det <- unlist(determinant(prec, logarithm = TRUE))["modulus"]
  return(- 0.5 * K * log(2 * pi) 
         - 0.5 * sum(log(diag(chol_prec_alpha))) 
         - 0.5 * p1 %*% t(p1) 
         - 0.5 * K * log(2 * pi) 
         - 0.5 * sum(log(diag(chol_prec_beta_1))) 
         - 0.5 * p2 %*% t(p2) 
         - 0.5 * K * log(2 * pi) 
         - 0.5 * sum(log(diag(chol_prec_beta_2))) 
         - 0.5 * p3 %*% t(p3)
         + sum(dgamma(sigma,
                      shape = shape_sigma,
                      rate = rate_sigma,
                      log = TRUE
         )))
}
# log_q(lambda, theta[1,])

# Algorithm ---------------------------------------------------------------

lambda_out <- run_VB(lambda)

# Check results
mean_alpha <- lambda_out[lambda_index$mu_alpha]
mean_beta_1 <- lambda_out[lambda_index$mu_beta[1,]]
mean_beta_2 <- lambda_out[lambda_index$mu_beta[2,]]
# The first alpha and beta should be a bit bigger

theta_out <- sim_theta(3600, lambda_out)

alpha <- colMeans(theta_out[,1:K])
beta_1 <- colMeans(theta_out[,(K + 1):(2 * K)])
beta_2 <- colMeans(theta_out[,(2 * K +1):(3 * K)])
sigma <- colMeans(theta_out[,(3 * K + 1):(3 * K + 2)])
f <- matrix(NA, ncol = K, nrow = n)
for (i in 1:n) {
  for (k in 1:K) {
    f[i, k] <- alpha[k] + beta_1[k] * x1[i] + beta_2[k] * x2[i]
  }
}
p <- matrix(NA, ncol = K, nrow = n)
for (i in 1:n) {
  p[i, ] <- exp(f[i, 1:K]) / (sum((exp(f[i, 1:K]))))
}
# 
# 
# f<-matrix(NA, ncol = K, nrow = n)
# for(i in 1:n){
#   for(k in 1:K){
#     f[i,k]<-mean_alpha[k] + mean_beta_1[k] * x1[i] + mean_beta_2[k] * x2[i]
#   }
# }
# p<-matrix(NA, ncol = K, nrow = n)
# for(i in 1:n){
#   p[i,] <- exp(f[i,1:K]) / (sum((exp(f[i,1:K]))))
# }
