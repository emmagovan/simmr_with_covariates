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
mu_s <- sources[, c(2, 3)] #+ disc[, c(2,4)]
sigma_s <- sources[, c(4, 5)]
mu_c <- TEFs[, c(2, 3)]
sigma_c <- TEFs[, c(4, 5)]
q <- conc[, c(2:3)]
#consumer$Skull, consumer$Wing, consumer$`Net Wt`
x <- matrix(c(rep(1,9),consumer$Sex, consumer$Age, consumer$Skull), 
            nrow = 4, 
            byrow = TRUE)
n_covariates <- nrow(x) -1
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

# Get the data into simmr
simmr_groups = simmr_load(mixtures=as.matrix(y),
                          source_names=unlist(sources[,1]),
                          source_means=as.matrix(sources[,c(2,3)]),
                          source_sds=as.matrix(sources[,c(4,5)]),
                          correction_means=as.matrix(TEFs[,c(2,3)]),
                          correction_sds=as.matrix(TEFs[,c(4,5)]),
                          concentration_means = as.matrix(conc[,2:3]))


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
c_0 <- c(rep(0.01, n_isotopes)) #Change to 0.0001
d_0 <- c(rep(0.01, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_covariates +1),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

x_scaled <- scale(x)


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
  beta <- matrix(theta[1:((n_covariates +1) * K)], nrow = (n_covariates+1), byrow = TRUE)
  sigma <- theta[((n_covariates +1) * K +1):(((n_covariates +1) * K)+n_isotopes)]
  f <- matrix(NA, ncol = K, nrow = n) 
  
  #Need to double check that this maths is right!!
  
   for (i in 1:n) {
     for (k in 1:K) {
      f[i,k] <-  sum(x_scaled[,i] * beta[,k])
     }
   }
  p <- matrix(NA, ncol = K, nrow = n)
  
  for (i in 1:n) {
    p[i, ] <- exp(f[i,]) / (sum((exp(f[i,]))))
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
  
  shape_sigma <- lambda[lambda_index$c]
  rate_sigma <- lambda[lambda_index$d]
  
  # Extract alpha, beta and sigma from theta
  beta <- matrix(theta[1:((n_covariates +1) * K)], nrow = (n_covariates+1), ncol = K,  byrow = TRUE)
  sigma <- theta[((n_covariates +1) * K +1):(((n_covariates +1) * K)+n_isotopes)]
  
  # This is a placeholder for more advanced code using chol_prec directly
  # prec <- crossprod(chol_prec)
  p_mat <- matrix(NA, nrow = n_covariates + 1, ncol = K) #row for each beta
  for(l in 1:(n_covariates +1)){
    p_mat[l,] <- (matrix(beta[l,] - mean_beta[l,], nrow = 1) %*% t(chol_prec_beta[,,l]))
  }
  
 sum_p = 0
 for(l in 1:(n_covariates +1)){
   sum_p = sum_p - 0.5 * K * log(2 * pi)-
                   0.5 * sum(log(diag(chol_prec_beta[,,l])))-
                   0.5 * matrix(p_mat[l,], nrow = 1) %*% (p_mat[l,])
                  
 }
  
  return(sum_p
         + sum(dgamma(sigma,
                      shape = shape_sigma,
                      rate = rate_sigma,
                      log = TRUE
         )))
}
# log_q(lambda, theta[1,])

# Algorithm ---------------------------------------------------------------

lambda_out <- run_VB(lambda)

n_samples <- 3600

# Check results ---------------------------------
theta_out <- sim_theta(n_samples, lambda_out)

#Easy way
beta<-matrix(colMeans(theta_out[,1:(K*(n_covariates+1))]), ncol = (n_covariates +1), byrow = TRUE)
sigma <- colMeans(theta_out[,(K*(n_covariates+1)+1):(K*(n_covariates+1)+n_isotopes)])

f1 <- matrix(NA, ncol = K, nrow = n)
 for (i in 1:n) {
   for (k in 1:K) {
     f1[i,k] <- x_scaled[1,i] * beta[1,k] + 
       x_scaled[2,i] * beta[2,k]
     x_scaled[3,i] * beta[3,k]
     x_scaled[4,i] * beta[4,k]
   
   }
 }


p1 <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1[i, ] <- exp(f1[i, 1:K]) / (sum((exp(f1[i, 1:K]))))
}





#### Making multiple samples


beta <- matrix(theta_out[,1:((n_covariates +1) * K)], ncol = (K *(n_covariates+1)), nrow = n_samples,  byrow = TRUE)
# beta <- matrix(theta_out[1,1:((n_covariates +1) * K)], nrow = (n_covariates+1), byrow = TRUE)
# sigma <- theta_out[1, ((n_covariates +1) * K +1):(((n_covariates +1) * K)+n_isotopes)]
#beta<-matrix(lambda_out[1:(n_covariates+1 * K)]

f <- array(NA, dim = c(n, K, n_samples))
#X<-(rbind(c(rep(1, n)), x_scaled)) #This just adds a row of ones for the alpha


for(s in 1:n_samples){
for (i in 1:n) {
  for (k in 1:K) {
    f[i,k,s] <-  sum(x_scaled[,i] * beta[s,k])
  }
}
}

p <- array(NA, dim = c(n, K, n_samples))

for (i in 1:n) {
  for (s in 1:n_samples){
  p[i,,s ] <- exp(f[i,,s]) / (sum((exp(f[i,,s]))))
  }
}










