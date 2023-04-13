library(tidyverse)
library(readxl)
library(R2jags)
library(simmr)

path = system.file("extdata", "geese_data.xls", package = "simmr")
geese_data = lapply(excel_sheets(path), read_excel, path = path)
sex <- geese_data[[1]]$Sex[1:9]
consumer <- read.csv("geese_consumer.csv", header = TRUE, stringsAsFactors = FALSE)
sources <- read.csv("geese_sources.csv", header = TRUE, stringsAsFactors = FALSE)[1:4,]
disc <- read.csv("geese_discrimination.csv", header = TRUE, stringsAsFactors = FALSE)
y<-(consumer %>% filter(Group == "1") %>% select(d15N, d13C))[1:9,]
c_0 <- c(1,1)
d_0 <- c(1,1)
n_isotopes <- 2

mu_kj <- sources[,c(2,4)] #+ disc[, c(2,4)]
K<-nrow(mu_kj)
sigma_kj <- sources[,c(3,5)]
mu_c <- disc[, c(2,4)]
sigma_c <- disc[, c(3,5)]
q <- sources[, c(6:7)]


simmr_obj <- simmr_load(y, 
                        sources$Sources,
                        mu_kj,
                        sigma_kj,
                        mu_c,
                        sigma_c,
                        q)

sim_out <-simmr_mcmc(simmr_obj)


summary(sim_out, type = "statistics")



