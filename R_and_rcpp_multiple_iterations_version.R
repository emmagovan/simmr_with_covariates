#R and RCPP code - can swap from one to another
#Multiple iterations version

#This just reads in the geese data, the 
source("data.R")

#Only things we might want to change
c_0 <- c(rep(0.1, n_isotopes)) #Change to 0.001
d_0 <- c(rep(0.1, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))

lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

source("r_fns.R")
Rcpp::sourceCpp("run_VB.cpp")

run_VB_cpp_r_hybrid <- function(lambda, # Starting value of lambda
                                S = 100, # Number of samples to take
                                P = 10, # Maximum patience before you stop
                                beta_1 = 0.9, # Learning rates
                                beta_2 = 0.9, # Learning rates
                                tau = 100, # Iteration at which learning rate starts to decrease
                                eps_0 = 0.1, # Raw learning rate multiplier
                                t_W = 50 # Time window for working out convergence
) {
  
  # Starting
  set.seed(12345)
  theta <- sim_theta(S, lambda)
  # theta <- sim_thetacpp(S, lambda, K, n_isotopes, n_covariates)
  
  c <- control_var(lambda, theta)
  #c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       # 0.001, n_covariates, x_scaled, concentrationmeans,
                       # sourcemeans, correctionmeans, corrsds,sourcesds,y)
  
  g_0 <- nabla_LB(lambda, theta)
  
  #g_0 <- nabla_LB_cpp(lambda,  theta, 
                      # K, n_isotopes, beta_prior,
                      # S,  n_covariates,
                      # x_scaled,
                      # concentrationmeans,  sourcemeans,
                      # correctionmeans,
                      # corrsds,  sourcesds,  y,
                      # rep(0, length(lambda)))
  nu_0 <- g_0^2
  g_bar <- g_0
  nu_bar <- nu_0
  
  # Set up
  t <- 1
  patience <- 0
  stop <- FALSE
  LB <- rep(NA, t_W)
  max_LB_bar <- -Inf
  
  while (!stop) {
    if (t %% 10 == 0) print(t)
    print("t")
    print(t)
    # Generate new samples
    set.seed(12344)
    theta <- sim_theta(S, lambda)
    #theta <-sim_thetacpp(S, lambda, K, n_isotopes, n_covariates)
    
    # Compute g_t
    
    g_t <- nabla_LB(lambda, theta, c)
    # g_t <- nabla_LB_cpp(lambda,  theta, 
    #                     K, n_isotopes, beta_prior,
    #                     S,  n_covariates,
    #                     x_scaled,
    #                     concentrationmeans,  sourcemeans,
    #                     correctionmeans,
    #                     corrsds,  sourcesds,  y,
    #                     c)
    
    print("g_t")
    print(g_t)
    
    # Compute new control variate
    
    c<- control_var(lambda, theta)
    # c <- control_var_cpp(lambda, theta, K, n_isotopes,
    #                      0.001, n_covariates, x_scaled, concentrationmeans,
    #                      sourcemeans, correctionmeans, corrsds,sourcesds,y)
    

    print("c")
    print(c)
    
    # Update the learning rates
    nu_t <- g_t^2
    print("nu_t")
    print(nu_t)
    g_bar <- beta_1 * g_bar + (1 - beta_1) * g_t
    nu_bar <- beta_2 * nu_bar + (1 - beta_2) * nu_t
    
    # Update the learning rate
    alpha_t <- min(eps_0, eps_0 * tau / t)
    print("alpha_t")
    print(alpha_t)
    
    # Update lambda
    lambda <- lambda + alpha_t * g_bar / sqrt(nu_bar)
    print("lambda")
    print(lambda)
    
    # Compute the moving average LB if out of warm-up
    if (t <= t_W) {
      # Compute a new lower bound estimate
      # LB[t] <- LB_lambda_cpp( theta,  lambda, 
      #                         hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
      #                         beta_prior,
      #                         n_covariates,
      #                         x_scaled,
      #                         concentrationmeans,  sourcemeans,
      #                         correctionmeans,
      #                         corrsds,  sourcesds,  y)
      LB[t] <- LB_lambda(lambda, theta)
    } else {
      LB[1:(t_W - 1)] <- LB[2:t_W]
      LB[t_W] <- LB_lambda(lambda, theta)
      # LB[t_W] <- LB_lambda_cpp(theta,  lambda, 
      #                          hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
      #                          beta_prior,
      #                          n_covariates,
      #                          x_scaled,
      #                          concentrationmeans,  sourcemeans,
      #                          correctionmeans,
      #                          corrsds,  sourcesds,  y)
      LB_bar <- mean(LB)
      max_LB_bar <- max(max_LB_bar, LB_bar)
      if (LB_bar >= max_LB_bar) {
        patience <- 0
      } else {
        patience <- patience + 1
      }
    }
    
    if (patience > P) {
      print("Completed!")
      stop <- TRUE
    }
    t <- t + 1
  }
  return(lambda)
}





lambdaout <- run_VB_cpp_r_hybrid(lambda)


n_samples = 3600

theta_out_rcpp <- sim_theta(n_samples, lambdaout)

beta_rcpp<-matrix(colMeans(theta_out_rcpp[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma_rcpp<- colMeans(theta_out_rcpp[,(n_sources*(n_covariates)+1):(n_sources*(n_covariates)+n_isotopes)])

f1_rcpp <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1_rcpp[,k] = x_scaled %*% beta_rcpp[,k]
}
p1_rcpp <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1_rcpp[i, ] <- exp(f1_rcpp[i, 1:K]) / (sum((exp(f1_rcpp[i, 1:K]))))
}

print(p1_rcpp)





run_VB_cpp_r_hybrid <- function(lambda, # Starting value of lambda
                                S = 100, # Number of samples to take
                                P = 10, # Maximum patience before you stop
                                beta_1 = 0.9, # Learning rates
                                beta_2 = 0.9, # Learning rates
                                tau = 100, # Iteration at which learning rate starts to decrease
                                eps_0 = 0.1, # Raw learning rate multiplier
                                t_W = 50 # Time window for working out convergence
) {
  
  # Starting
  set.seed(12345)
  theta <- sim_theta(S, lambda)
  #theta <- sim_thetacpp
  c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       0.001, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds,sourcesds,y)
  g_0 <- nabla_LB_cpp(lambda,  theta, 
                      K, n_isotopes, beta_prior,
                      S,  n_covariates,
                      x_scaled,
                      concentrationmeans,  sourcemeans,
                      correctionmeans,
                      corrsds,  sourcesds,  y,
                      rep(0, length(lambda)))
  nu_0 <- g_0^2
  g_bar <- g_0
  nu_bar <- nu_0
  
  # Set up
  t <- 1
  patience <- 0
  stop <- FALSE
  LB <- rep(NA, t_W)
  max_LB_bar <- -Inf
  
  while (!stop) {
    if (t %% 10 == 0) print(t)
    print("t")
    print(t)
    # Generate new samples
    set.seed(12345)
    theta <- sim_theta(S, lambda)
    
    # Compute g_t
    g_t <- nabla_LB_cpp(lambda,  theta, 
                        K, n_isotopes, beta_prior,
                        S,  n_covariates,
                        x_scaled,
                        concentrationmeans,  sourcemeans,
                        correctionmeans,
                        corrsds,  sourcesds,  y,
                        c)
    print("g_t")
    print(g_t)
    
    # Compute new control variate
    c <- control_var_cpp(lambda, theta, K, n_isotopes,
                         0.001, n_covariates, x_scaled, concentrationmeans,
                         sourcemeans, correctionmeans, corrsds,sourcesds,y)
    
    # c <- control_var(lambda, theta)
    print("c")
    print(c)
    
    # Update the learning rates
    nu_t <- g_t^2
    print("nu_t")
    print(nu_t)
    g_bar <- beta_1 * g_bar + (1 - beta_1) * g_t
    nu_bar <- beta_2 * nu_bar + (1 - beta_2) * nu_t
    
    # Update the learning rate
    alpha_t <- min(eps_0, eps_0 * tau / t)
    print("alpha_t")
    print(alpha_t)
    
    # Update lambda
    lambda <- lambda + alpha_t * g_bar / sqrt(nu_bar)
    print("lambda")
    print(lambda)
    
    # Compute the moving average LB if out of warm-up
    if (t <= t_W) {
      # Compute a new lower bound estimate
      LB[t] <- LB_lambda_cpp( theta,  lambda, 
                              hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                              beta_prior,
                              n_covariates,
                              x_scaled,
                              concentrationmeans,  sourcemeans,
                              correctionmeans,
                              corrsds,  sourcesds,  y)
    } else {
      LB[1:(t_W - 1)] <- LB[2:t_W]
      LB[t_W] <- LB_lambda_cpp(theta,  lambda, 
                               hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                               beta_prior,
                               n_covariates,
                               x_scaled,
                               concentrationmeans,  sourcemeans,
                               correctionmeans,
                               corrsds,  sourcesds,  y)
      LB_bar <- mean(LB)
      max_LB_bar <- max(max_LB_bar, LB_bar)
      if (LB_bar >= max_LB_bar) {
        patience <- 0
      } else {
        patience <- patience + 1
      }
    }
    
    if (patience > P) {
      print("Completed!")
      stop <- TRUE
    }
    t <- t + 1
  }
  return(lambda)
}



lambdaout <- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.1, 
                             as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                             as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                             100, 10, 0.9, 0.9, 1000, 0.1, 50)


n_samples = 3600

theta_out_rcpp <- sim_thetacpp(n_samples, lambdaout)

beta_rcpp<-matrix(colMeans(theta_out_rcpp[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma_rcpp<- colMeans(theta_out_rcpp[,(n_sources*(n_covariates)+1):(n_sources*(n_covariates)+n_isotopes)])

f1_rcpp <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1_rcpp[,k] = x_scaled %*% beta_rcpp[,k]
}
p1_rcpp <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1_rcpp[i, ] <- exp(f1_rcpp[i, 1:K]) / (sum((exp(f1_rcpp[i, 1:K]))))
}

print(p1_rcpp)





