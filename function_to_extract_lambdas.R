# function to extract lambdas


lambda_index <- function(n_covariates, K, n_isotopes){
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
