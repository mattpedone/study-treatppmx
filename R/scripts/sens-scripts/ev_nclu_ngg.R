rm(list = ls())
library(vweights)
library(treatppmx)
set.seed(121)

nobs <- 75
#beta <- 10
sigma <- .0
kappa <- 1

w <- c(round(vweights::computepnclu(nobs, sigma, kappa), 5))
v <- 1:length(w)
ev <- sum(w*v)
ev