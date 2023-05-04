#JAGS
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
x <- matrix(c(rep(1,9), consumer$Sex, consumer$Age, consumer$Skull), 
            nrow = 4, 
            byrow = TRUE)

# Xmat<-matrix(c(rep(1,9),consumer$Sex, consumer$Wing, consumer$Skull, consumer$`Net Wt`), 
#              nrow = 5, 
#              byrow = TRUE)
Xmat<-scale(x)

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
      mu_f[i,k] = inprod(x[,i], beta[,k])
  }
  
  
  for(a in 1:(ncov +1)){
   beta[a,k] ~ dnorm(0,1)
  }
  }   
  
  for(j in 1:J) { 
      sigma[j] ~ dgamma(1, 1)
      }
}
"

model_data = list(y=consumer,
                  s_mean=mu_s,
                  s_sd=sigma_s,
                  c_mean=mu_c,
                  c_sd=sigma_c,
                  q=q, 
                  N=nrow(consumer),
                  K=4,
                  ncov = nrow(Xmat) -1,
                  J=2,
                  x = Xmat)


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



