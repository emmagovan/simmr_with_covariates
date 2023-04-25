#general number of covariates
#Specify the number at the start and then automatically generate it??

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
x <- matrix(c(consumer$Sex, consumer$Wing, consumer$Skull, consumer$`Net Wt`), 
            nrow = 4, 
            byrow = TRUE)
n_covariates <- nrow(x)
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

# Variables for the FFVB
S <- 100
mu_beta_zero <- matrix(c(rep(0, K * (n_covariates +1))), 
                       nrow = (n_covariates +1), 
                       ncol = K)
#n_covariates + 1 for alpha (or beta_0) and then 1 beta for each covariate
sigma_beta_zero <- matrix(c(rep(1, K * (n_covariates +1))), 
                          nrow = (n_covariates +1), 
                          ncol = K)

n_isotopes <- ncol(mu_c)
c_0 <- c(rep(1, n_isotopes))
d_0 <- c(rep(1, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_covariates +1),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

# function to extract lambdas --------------------------------------------
lambda_extract <- function(n_covariates, K, n_isotopes){
  mat_size = K * (K+1) /2
  mu_beta = matrix(data = NA, nrow = (n_covariates+1), ncol = K)
  sigma_beta = matrix(data = NA, nrow = (n_covariates+1), ncol = mat_size)
  
  for(i in 1:(n_covariates+1)){
    mu_beta[i,] = ((i-1) * mat_size + (i-1) * K +1):((i-1) * mat_size + (i-1) * K + K)
    sigma_beta[i,] = ((i-1) * mat_size + (i) * K +1): ((i-1) * mat_size + (i) * K +mat_size)
  }
  
  c = (sigma_beta[n_covariates+1, mat_size] + 1):(sigma_beta[n_covariates+1, mat_size] + n_isotopes)
  d = (sigma_beta[n_covariates+1, mat_size] + n_isotopes + 1):(sigma_beta[n_covariates+1, mat_size] + 2 * n_isotopes)
  
  return(list(mu_beta = mu_beta,
              sigma_beta = sigma_beta,
              c = c,
              d = d
  ))
}
#just use here for now - then build in to r
lambda_index <- lambda_extract(n_covariates, K, n_isotopes)


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
  
  # mean_alpha <- lambda[lambda_index$mu_alpha]
  # 
  # # K*(K-1) precision terms
  # chol_prec_alpha <- matrix(0, nrow = K, ncol = K)
  # chol_prec_alpha[upper.tri(chol_prec_alpha, diag = TRUE)] <-
  #   lambda[lambda_index$sigma_alpha]
  # 
  # mean_beta_1 <- lambda[lambda_index$mu_beta[1,]]
  # chol_prec_beta_1 <- matrix(0, nrow = K, ncol = K)
  # chol_prec_beta_1[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
  #   lambda[lambda_index$sigma_beta[1,]]
  # 
  # mean_beta_2 <- lambda[lambda_index$mu_beta[2,]]
  # chol_prec_beta_2 <- matrix(0, nrow = K, ncol = K)
  # chol_prec_beta_2[upper.tri(chol_prec_beta_1, diag = TRUE)] <-
  #   lambda[lambda_index$sigma_beta[2,]]
  # 
  ## Create a loop instead to do this I think? Will need to make an array?
  mean_beta <- matrix(0, nrow = (n_covariates +1), ncol = K)
  for(i in 1:(n_covariates+1)){
    mean_beta[i,] <- lambda[lambda_index$mu_beta[i,]]
  }
  
  chol_prec_beta <- array(data = 0, dim = c(K, K, (n_covariates +1)))
  
  for(i in 1:(n_covariates +1)){
    chol_prec_beta[,,i][upper.tri(chol_prec_beta[,,i], diag = TRUE)] <-
      lambda[lambda_index$sigma_beta[i,]]
  }
  a<-array(NA, dim =c(S, K, (n_covariates+1)))
  thetabeta<-matrix(NA, ncol = (n_covariates+1) * K, nrow = S)
  
  
  
  for(i in 1:(n_covariates+1)){
    
    thetabeta[,(1+(i-1)*K):((i)*K)] = t(rMVNormC(S, mu = mean_beta[i,], U = chol_prec_beta[,,i]))

  }
  
 

  
  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    thetabeta,
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
  # Create betas and sigma
  beta <- matrix(theta[1:((n_covariates +1) * K)], nrow = (n_covariates+1))
  sigma <- theta[((n_covariates +1) * K +1):(((n_covariates +1) * K)+n_isotopes)]
  f <- matrix(NA, ncol = K, nrow = n)
  X<-rbind(c(rep(1, n)), x) #This just adds a row of ones for the alpha
  
  #Need to double check that this maths is right!!
  
  for (i in 1:n) {
    for (k in 1:K) {
      f[i, ] <-  beta[k,] * X[k,i]
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
  
  beta_sum <- 0
  for (i in 1:(n_covariates+1)){
    beta_sum = beta_sum +sum(dnorm(beta[i,], mu_beta_zero[i,], sigma_beta_zero[i,], log = TRUE))
  }
  
  return(hold + beta_sum +
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



f<-matrix(NA, ncol = K, nrow = n)
for(i in 1:n){
  for(k in 1:K){
    f[i,k]<-mean_alpha[k] + mean_beta_1[k] * x1[i] + mean_beta_2[k] * x2[i]
  }
}
p<-matrix(NA, ncol = K, nrow = n)
for(i in 1:n){
  p[i,] <- exp(f[i,1:K]) / (sum((exp(f[i,1:K]))))
}
