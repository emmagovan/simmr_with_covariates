#For testing rcpp

library(Rcpp)

S = 100
n_sources = 4
n_tracers = 2
n_cov = 6

c_0 <- c(rep(1, n_tracers)) #Change to 0.0001
d_0 <- c(rep(1, n_tracers))
beta_lambda<-c(rep(0, n_sources),rep(17, n_sources * (n_sources + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_cov),
  rep(1, n_tracers), #shape
  rep(1, n_tracers) #rate
)


set.seed(123)
a<-sim_thetacpp(S, lambda, n_sources, n_tracers, n_cov)

n_covariates = n_cov -1
K = n_sources
b<-sim_theta(S, lambda)


b<-matrix(NA, nrow = 3, ncol = 4)

for(i in 1:3){
  for(j in 1:4){
    b[i,j] = (i-1)*4 + j
  }
}

n_cov = 2
n = 5
K = 4

theta<-c(1:(n_cov *K))
beta<-matrix(c(1:8), nrow = n_cov, ncol =K, byrow = TRUE)
x_scaled<-matrix(rnorm(n_cov*n,1,1), ncol = n_cov, nrow = n)

f<-x_scaled %*% beta



hfn(theta, K, n, n_cov, x_scaled)
