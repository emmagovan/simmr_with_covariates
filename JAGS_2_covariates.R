#JAGS

model_code <- "
model{
  for (i in 1:N) {
    for (j in 1:J) { 
      y[i,j] ~ dnorm(inprod(p[i,]*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p[i,],q[,j]), 1/var_y[i,j])
      var_y[i,j] <- inprod(pow(p[i,]*q[,j],2),s_sd[,j]^2+c_sd[,j]^2)/pow(inprod(p[i,],q[,j]),2) 
        + pow(sigma[j],2)
    }
  }
  for(i in 1:N) {
    p[i,1:K] <- expf[i,]/sum(expf[i,])
    for(k in 1:K) {
      expf[i,k] <- exp(f[i,k])
      f[i,k] = mu_f[i,k]
    }
  }
  for(k in 1:K) { 
  for(i in 1:N) {
      mu_f[i,k] <- alpha[k] +beta1[k]*x1[i] +beta2[k]*x2[i]
      }
      
    alpha[k] ~ dnorm(0,1)
    beta1[k] ~dnorm(0,1) 
    beta2[k] ~dnorm(0,1) 
  }

  
  for(j in 1:J) { sigma[j] ~ dgamma(0.001, 0.001) }
}
"

model_data = list(y=y,
                  s_mean=mu_kj,
                  s_sd=sigma_kj,
                  c_mean=mu_c,
                  c_sd=sigma_c,
                  q=q, 
                  N=nrow(y),
                  K=length(sources$Sources),
                  J=ncol(y),
                  x1 = x1,
                  x2 = x2)


model_parameters <- c("alpha", "beta", "sigma", "p")

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
```