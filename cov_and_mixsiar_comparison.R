#test run_VB.cpp with multiple covariates, different data sets etc

# Load in package
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

# Source in the generic functions
Rcpp::sourceCpp("run_VB.cpp")


# Extract data ------------------------------------------------------------

path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] #|> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]

# Put them into the right names
n <- nrow(consumer)
n_isotopes <- 2
K <- nrow(sources)
mu_s <- as.matrix(sources[, c(2, 3)] )
sigma_s <- as.matrix(sources[, c(4, 5)])
mu_c <- as.matrix(TEFs[, c(2, 3)])
sigma_c <- as.matrix(TEFs[, c(4, 5)])
q <- as.matrix(conc[, c(2:3)])



########## SET THESE
groups <- simmr::geese_data$groups
groups_vec = c(rep(1,9), rep(2,29), rep(3,74), rep(4,10), 
               rep(5,41), rep(6,20), rep(7,32), rep(8,36))

x <- matrix(c((groups_vec)),
ncol = 1)
x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))
#x_scaled <- matrix(c(rep(1, n)), ncol = 1)

n_covariates <- (ncol(x_scaled))
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

#-------------------- FFVB variables-----------------------
#Don't actually need to edit anything here, just run it


# Variables for the FFVB
S <- 100
mu_beta_zero <- matrix(c(rep(0, K * (n_covariates))), 
                       nrow = (n_covariates), 
                       ncol = K)
#n_covariates + 1 for alpha (or beta_0) and then 1 beta for each covariate
sigma_beta_zero <- matrix(c(rep(1, K * (n_covariates))), 
                          nrow = (n_covariates), 
                          ncol = K)

n_isotopes <- ncol(mu_c)
c_0 <- c(rep(0.1, n_isotopes)) #Change to 0.0001
d_0 <- c(rep(0.1, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)


lambdaout <- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.1, 
                        as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                        as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                        100, 10, 0.9, 0.9, 100, 0.1, 50)


n_samples = 3600

theta_out_rcpp <- sim_thetacpp(S, lambdaout, K, n_isotopes, n_covariates)

beta_rcpp<-matrix(colMeans(theta_out_rcpp[,1:(K*(n_covariates))]), nrow = (n_covariates))
tau_rcpp<- colMeans(theta_out_rcpp[,(K*(n_covariates)+1):(K*(n_covariates)+n_isotopes)])

f1_rcpp <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1_rcpp[,k] = x_scaled %*% beta_rcpp[,k]
}
p1_rcpp <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1_rcpp[i, ] <- exp(f1_rcpp[i, 1:K]) / (sum((exp(f1_rcpp[i, 1:K]))))
}

print(p1_rcpp)



#--------------------JAGS--------------

model_code <- "
model{
  for (i in 1:N) {
    for (j in 1:J) { 
      y[i,j] ~ dnorm(inprod(p[i,]*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p[i,],q[,j]), 1/var_y[i,j])
      var_y[i,j] = inprod(pow(p[i,]*q[,j],2),s_sd[,j]^2+c_sd[,j]^2)/pow(inprod(p[i,],q[,j]),2) 
        + pow(sigma[j],2)
    }
  }
  for(i in 1:N) {
    p[i,1:K] = expf[i,]/sum(expf[i,])
    for(k in 1:K) {
      expf[i,k] = exp(f[i,k])
      f[i,k] = mu_f[i,k]
    }
  }
  for(k in 1:K) {
  for(i in 1:N) {
      mu_f[i,k] = inprod(x[i,], beta[,k])
  }
  
  
  for(a in 1:(ncov)){
   beta[a,k] ~ dnorm(0,1)
  }
  }   
  
  for(j in 1:J) { 
      sigma[j] ~ dgamma(0.1, 0.1)
      }
}
"

model_data = list(y=y,
                  s_mean=as.matrix(mu_s),
                  s_sd=as.matrix(sigma_s),
                  c_mean=as.matrix(mu_c),
                  c_sd=as.matrix(sigma_c),
                  q=as.matrix(q), 
                  N=nrow(y),
                  K=4,
                  ncov = n_covariates,
                  J=2,
                  x = x_scaled)


model_parameters <- c("beta", "sigma", "p")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.chains = 4, # Number of different starting positions
  n.iter = 10000, # Number of iterations
  n.burnin = 2000, # Number of iterations to remove at start
  n.thin = 5
) # Amount of thinning)


print(model_run)



library(microbenchmark)
microbenchmark(run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 1, 
                          as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                          as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                          100, 10, 0.9, 0.9, 100, 0.1, 50),
               run_model(run="short", mix, source, discr, model_filename,
                         alpha.prior = 1, resid_err, process_err),
               times = 5L
               )




