# A simple version of the covariate model adapted from
# 'simplest_covariate_model.Rmd'

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
mu_s <- sources[, c(2, 4)] #+ disc[, c(2,4)]
sigma_s <- sources[, c(3, 5)]
mu_c <- TEFs[, c(2, 4)]
sigma_c <- TEFs[, c(3, 5)]
q <- conc[, c(2:3)]
x <- consumer$Sex
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

# Variables for the FFVB
S <- 100
mu_alpha_zero <- c(rep(0, K))
mu_beta_zero <- c(rep(0, K))
sigma_alpha_zero <- c(rep(1, K))
sigma_beta_zero <- c(rep(1, K))
n_isotopes <- ncol(mu_c)
c_0 <- c(rep(1, n_isotopes))
d_0 <- c(rep(1, n_isotopes))
lambda <- c(
  rep(0, K),
  rep(1, K * (K + 1) / 2),
  rep(0, K),
  rep(1, K * (K + 1) / 2),
  rep(1, n_isotopes),
  rep(1, n_isotopes)
)

# sim_theta ---------------------------------------------------------------

sim_theta <- function(S = 100, lambda) {
  # lambda contains the parameters
  # mean_alpha, 1:K
  # chol(Sigma_alpha), (K + 1):(K + (K * (K + 1)) / 2)
  # mean_beta, (K + (K * (K + 1) / 2) + 1):(2 * K + (K * (K + 1) / 2))
  # chol(Sigma_beta), (2 * K + (K * (K + 1) / 2) + 1):(K * (K + 1) + 2 * K)
  # sigma_j_shape, ((K * (K + 1) + 2 * K) + 1):((K * (K + 1) + 2 * K) + n_isotopes)
  # sigma_j_scale, ((K * (K + 1) + 2 * K) + n_isotopes + 1):((K * (K + 1) + 2 * K) + 2 * n_isotopes)

  mean_alpha <- lambda[1:K]

  # K*(K-1) precision terms
  chol_prec_alpha <- matrix(0, nrow = K, ncol = K)
  chol_prec_alpha[upper.tri(chol_prec_alpha, diag = TRUE)] <-
    lambda[(K + 1):(K + (K * (K + 1)) / 2)]

  mean_beta <- lambda[(K + (K * (K + 1) / 2) + 1):(2 * K + (K * (K + 1) / 2))]
  chol_prec_beta <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta[upper.tri(chol_prec_beta, diag = TRUE)] <-
    lambda[(2 * K + (K * (K + 1) / 2) + 1):(K * (K + 1) + 2 * K)]

  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    t(rMVNormC(S, mu = mean_alpha, U = chol_prec_alpha)),
    t(rMVNormC(S, mu = mean_beta, U = chol_prec_beta)), 
    matrix(
      rgamma(S * n_isotopes,
        shape = lambda[((K * (K + 1) + 2 * K) + 1):((K * (K + 1) + 2 * K) + n_isotopes)],
        rate = lambda[((K * (K + 1) + 2 * K) + n_isotopes + 1):((K * (K + 1) + 2 * K) + 2 * n_isotopes)]
      ),
      nrow = S,
      ncol = n_isotopes,
      byrow = TRUE
    )
  )

  return(theta)
}
theta <- sim_theta(S, lambda)
# Theta is alpha (K of these), beta (K of these), and sigma (J of these)
# lambda is mu_alpha, chol(sigma_alpha), mu_beta, chol(sigma_beta), c, d

# h -----------------------------------------------------------------------

# Log of likelihood added to prior
h <- function(theta) {
  # Create alpha, beta and sigma
  alpha <- theta[1:K]
  beta <- theta[(K + 1):(2 * K)]
  sigma <- theta[(2 * K + 1):(2 * K + 2)]
  f <- matrix(NA, ncol = K, nrow = n)
  for (i in 1:n) {
    for (k in 1:K) {
      f[i, k] <- alpha[k] + beta[k] * x[i]
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
    sum(dnorm(beta, mu_beta_zero, sigma_beta_zero, log = TRUE)) +
    sum(dgamma(sigma, shape = c_0, rate = d_0, log = TRUE)))
}
h(theta[1, ])


# log_q -------------------------------------------------------------------

log_q <- function(lambda, theta) {
  mean_alpha <- lambda[1:K]
  chol_prec_alpha <- matrix(0, nrow = K, ncol = K)
  chol_prec_alpha[upper.tri(chol_prec_alpha, diag = TRUE)] <-
    lambda[(K + 1):(K + (K * (K + 1)) / 2)]

  mean_beta <- lambda[(K + (K * (K + 1) / 2) + 1):(2 * K + (K * (K + 1) / 2))]
  chol_prec_beta <- matrix(0, nrow = K, ncol = K)
  chol_prec_beta[upper.tri(chol_prec_beta, diag = TRUE)] <-
    lambda[(2 * K + (K * (K + 1) / 2) + 1):(K * (K + 1) + 2 * K)]

  shape_sigma <- lambda[((K * (K + 1) + 2 * K) + 1):
                          ((K * (K + 1) + 2 * K) + n_isotopes)]
  rate_sigma <- lambda[((K * (K + 1) + 2 * K) + n_isotopes + 1):
                         ((K * (K + 1) + 2 * K) + 2 * n_isotopes)]
  
  # Extract alpha, beta and sigma from theta
  alpha <- theta[1:K]
  beta <- theta[(K + 1):(2 * K)]
  sigma <- theta[(2 * K + 1):(2 * K + 2)]
  
  # This is a placeholder for more advanced code using chol_prec directly
  # prec <- crossprod(chol_prec)
  p1 <- matrix(alpha - mean_alpha, nrow = 1) %*% t(chol_prec_alpha)
  p2 <- matrix(beta - mean_beta, nrow = 1) %*% t(chol_prec_beta)
  # log_det <- unlist(determinant(prec, logarithm = TRUE))["modulus"]
  return(- 0.5 * K * log(2 * pi) 
         - 0.5 * sum(log(diag(chol_prec_alpha))) 
         - 0.5 * p1 %*% t(p1) 
         - 0.5 * K * log(2 * pi) 
         - 0.5 * sum(log(diag(chol_prec_beta))) 
         - 0.5 * p2 %*% t(p2) 
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
mean_alpha <- lambda[1:K]
mean_beta <- lambda[(K + (K * (K + 1) / 2) + 1):(2 * K + (K * (K + 1) / 2))]
# The first alpha and beta should be a bit bigger