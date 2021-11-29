rm(list = ls())
library(vweights)
library(treatppmx)
set.seed(121)

nobs <- 75 
beta <- 1
sigma <- .01
kappa <- beta*sigma; kappa

w <- c(round(vweights::computepnclu(nobs, sigma, kappa), 5))
v <- 1:length(w)
ev <- sum(w*v)
ev