rm(list = ls())
library(vweights)
library(treatppmx)
set.seed(121)

nobs <- 50 
beta <- 48
sigma <- .25
kappa <- beta*sigma; kappa

w <- c(round(vweights::computepnclu(nobs, sigma, kappa), 5))
v <- 1:length(w)
ev <- sum(w*v)
ev