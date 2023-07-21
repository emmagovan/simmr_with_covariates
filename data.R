#Data file
#Rcpp code running one iteration at a time and printing everything
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

# Extract data ------------------------------------------------------------
path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(            excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] |> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]

# Put them into the right names
n <- nrow(consumer)
n_isotopes <- 2
K <- nrow(sources)
n_sources = K
mu_s <- sources[, c(2, 3)] 
sigma_s <- sources[, c(4, 5)]
mu_c <- TEFs[, c(2, 3)]
sigma_c <- TEFs[, c(4, 5)]
q <- conc[, c(2:3)]

concentrationmeans = as.matrix(q)
sourcemeans = as.matrix(mu_s)
correctionmeans = as.matrix(mu_c)
corrsds = as.matrix(sigma_c)
sourcesds = as.matrix(sigma_s)
beta_prior = 0.1

########## SET THESE
#   x <- matrix(c(consumer$Skull), 
# ncol = 1)
# x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))
x_scaled <- matrix(c(rep(1, 9)), ncol = 1)

n_covariates <- (ncol(x_scaled))
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()
# y <- consumer |>
#   select(d13C_Pl) |>
#   as.matrix()

# Get the data into simmr
simmr_groups = simmr_load(mixtures=as.matrix(y),
                          source_names=unlist(sources[,1]),
                          source_means=as.matrix(sources[,c(2,3)]),
                          source_sds=as.matrix(sources[,c(4,5)]),
                          correction_means=as.matrix(TEFs[,c(2,3)]),
                          correction_sds=as.matrix(TEFs[,c(4,5)]),
                          concentration_means = as.matrix(conc[,2:3]))


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
c_0 <- c(rep(0.1, n_isotopes)) #Change to 0.001
d_0 <- c(rep(0.1, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))

lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

