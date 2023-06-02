#For testing rcpp

library(Rcpp)
library(simmr)
library(tidyverse)

path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] |> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]


concmeans<-conc[, c(2:3)]
corrmeans<-TEFs[, c(2, 3)]
sourcemeans<-sources[, c(2, 3)]
sourcesds<-sources[, c(4, 5)]
corrsds <- TEFs[, c(4, 5)]
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

S = 100
n_sources = 4
n_tracers = 2
n_cov = 6

c_0 <- c(rep(1, n_tracers)) #Change to 0.0001
d_0 <- c(rep(1, n_tracers))
beta_lambda<-c(rep(0, n_sources),rep(1, n_sources * (n_sources + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_cov),
  rep(1, n_tracers), #shape
  rep(1, n_tracers) #rate
)

x <- matrix(c(consumer$Sex, consumer$Skull, consumer$Wing, consumer$Age, (consumer$`Net Wt`/100)), 
            ncol = 5)
x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))

set.seed(123)
theta<-sim_thetacpp(S, lambda, n_sources, n_tracers, n_cov)

n_covariates = n_cov -1
K = n_sources
b<-sim_theta(S, lambda)


b<-matrix(NA, nrow = 3, ncol = 4)

for(i in 1:3){
  for(j in 1:4){
    b[i,j] = (i-1)*4 + j
  }
}




p<-hfn(theta, K, n, n_cov, x_scaled)

theta1 <- c(0.1659410, -1.2335401, -1.3089230, 2.4150564, -0.3424654, 0.2279773,
             0.1665621, -0.2158007, -0.7411550, -1.1588777,  2.2865975, -1.0242991,
              0.8267962,  1.1047661, -1.5606877,  0.7538129,  0.6451374, -0.1970286,
             0.5416076,  0.1610883, -0.7182411,  3.2600013, -2.5644053,  1.2555209,
             0.9267664,  0.4518883)


a<-hcpp(4, 2, 6, 0.001, x_scaled, 
     as.matrix(concmeans), as.matrix(sourcemeans), as.matrix(corrmeans), 
     as.matrix(corrsds), as.matrix(sourcesds), theta1, y)

#comparison
for(j in 1:n_tracers){
for(i in 1:n){
    print(sqrt(sum(p[i, ]^2 * concmeans[, j]^2 * (sourcemeans[, j]^2 + corrmeans[, j]^2))/sum(p[i,]^2*concmeans[,j]^2)))
  }
}


log_q_cpp(theta1, lambda, 4, 2, 100, 6)
       








