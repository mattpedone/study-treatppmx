rm(list=ls())
#set.seed(121)
load("data/LGGdata.rda")
library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

#name <- c("a1s01.RData")
trtsgn <- c(matchRTComp[,10]) + 1
npat <- length(trtsgn)

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- rep(0, 6)
gof_all <- rep(0, 2)
#myres0 <- sellines_all <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)

cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

X <- data.frame(scale(matchRTComp[,16:38]))
Z <- data.frame(scale(matchRTComp[,c(11,13)]))#data.frame(orgx)#
#aggiusta!
Y <- matrix(0, nrow = nrow(X), ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))

n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 15000#0; 
burnin <- 7500#0; 
thinning <- 10

nout <- (iterations-burnin)/thinning
predAPT <- c()

sub <- sample(1:158, 1)

out_ppmx <- ppmxct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]), 
                   Xpred = data.frame(X[sub,]), Z = data.frame(Z[-sub,]), 
                   Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                   PPMx = 1, cohesion = 2, alpha = 1, sigma = 0.25,
                   similarity = 2, consim = 1, similparam = vec_par, 
                   calibration = 1, coardegree = 2, modelpriors, update_hierarchy = T,
                   hsp = T, iter = iterations, burn = burnin, thin = thinning, 
                   mhtunepar = c(0.05, 0.05), CC = n_aux, reuse = 1, nclu_init = 5)
      
#number of a cluster, mean, binder &varinf ----
mc <- apply(out_ppmx$nclu, 1, mean)
trt <- trtsgn[-sub]
num_treat <- table(trt)

cls1 <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
psm1 <- comp.psm(cls1)
mc_b1 <- minbinder.ext(psm1); 
mc_vi1 <- minVI(psm1); 

cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
psm2 <- comp.psm(cls2)
mc_b2 <- minbinder.ext(psm2); 
mc_vi2 <- minVI(psm2); 

mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))

mc_b; mc_vi

#posterior predictive probabilities ----
A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml);#A0

  #sellines <- 1:npat
A1 <- A0[1:3]%*%wk 
A2 <- A0[4:6]%*%wk 

A0; A1-A2
