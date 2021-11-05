library(parallel)
library(doParallel)
library(dplyr)
library(tibble)
library(tidyselect)
library(ggplot2)
library(vweights)
library(treatppmx)
set.seed(11)

K <- 10
#X <- data.frame(simupats)
#X <- X[,c(1:3)]

gendata <- function(n = 1000, pro = c(0.2,0.5,0.3), dim = 3){
  if(length(pro) != 3){
    stop("The length of probabilities must be equal to 3!")
  }
  data <- matrix(nrow = n, ncol = dim)
  for(i in 1:n){
    u <- runif(1)
    if(u < pro[1]){
      #cat(i,"ciao1","\n")
      data[i,] <- mvtnorm::rmvnorm(1, mean = rep(-2.1, dim), sigma = diag(.5, dim))
    }
    else{
      if(u < (pro[1] + pro[2])){
        #cat(i,"ciao2","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(0, dim), sigma = diag(.5, dim))
      }
      else{
        #cat(i,"ciao3","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(2.3, dim), sigma = diag(.5, dim))
      }
    }
  }
  return (data)
}

X <- gendata(n=100, dim = 10)

nobs <- nrow(X)

vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 10000#500
burnin <- 2500#00
thinning <- 10

nout <- (iterations-burnin)/thinning

cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

res_121 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = 2,
                                 sigma = .1, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 0,
                                 coardegree = 1, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_122 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = 2,
                                 sigma = .1, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 1,
                                 coardegree = 2, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_221 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = 2,
                                 sigma = .1, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 2,
                                 coardegree = 1, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_222 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = 2,
                                 sigma = .1, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 2,
                                 coardegree = 2, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_0 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 0, cohesion = 1, alpha = 2,
                                 sigma = .1, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 0,
                                 coardegree = 1, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res <- t(rbind(apply(res_121/nout, 2, mean), apply(res_122/nout, 2, mean), 
               apply(res_221/nout, 2, mean), apply(res_222/nout, 2, mean), 
               apply(res_0/nout, 2, mean), apply(res_121/nout, 2, sd), 
               apply(res_122/nout, 2, sd), apply(res_221/nout, 2, sd), 
               apply(res_222/nout, 2, sd), apply(res_0/nout, 2, sd), 
               c(round(vweights::computepnclu(100, .1, 2), 5))))

res <- (res[which(rowSums(res) != 0),])
res <- rownames_to_column(as.data.frame(res), var = "cluster")
newres <- cbind(similarity = c(rep("nocal", nrow(res)), rep("cal", nrow(res)), 
                               rep("coa1", nrow(res)), rep("coa2", nrow(res)), 
                               rep("ppm", nrow(res)), rep("ngg", nrow(res))),
                cluster = rep(res[,1], 6),
                freq = c(unlist(res[, c(2:6, 12)])),
                sd = c(unlist(res[, c(7:11)]), rep(0, nrow(res))))

dfres <- newres %>%
  as_tibble() %>% # lo trasformo in un dataframe
  #tibble::rownames_to_column(var = "cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  mutate(freq = as.numeric(freq)) %>%
  mutate(sd = as.numeric(sd))



ggplot(dfres, aes(x=cluster, y=freq, group=similarity, color=similarity)) +
  geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.1) +
  geom_line() +
  geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_minimal()

### ----
#gendata <- function(n = 1000, pro = c(0.2,0.5,0.3)){
#  if(length(pro) != 3){
#    stop("The length of probabilities must be equal to 3!")
#  }
#  data <- matrix(nrow = n, ncol = 3)
#  for(i in 1:n){
#    u <- runif(1)
#    if(u < pro[1]){
#      #cat(i,"ciao1","\n")
#      data[i,] <- rmvnorm(1,mean=rep(-2.1, 3) ,sigma=diag(0.5, 3))
#    }
#    else{
#      if(u < (pro[1] + pro[2])){
#        #cat(i,"ciao2","\n")
#        data[i,] <- rmvnorm(1,mean=rep(0, 3),sigma=diag(0.5, 3))
#      }
#      else{
#        #cat(i,"ciao3","\n")
#        data[i,] <- rmvnorm(1,mean=rep(2.3, 3),sigma=diag(0.5, 3))
#      }
#    }
#  }
#  return (data)
#}
#
#X <- gendata(n=100)

