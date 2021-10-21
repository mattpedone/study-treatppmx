rm(list=ls())
load("data/modscenario2.rda")
library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

set.seed(121)

name <- c("NGG_alpha_1sigma5.RData")
K <- k <- 1 #repliche
npat <- length(trtsgn)

predAPT_all <- array(0, dim = c(npat, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)
sellines_all <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)


cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

#X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
X <- data.frame(mydata)
Z <- data.frame(cbind(myz2, myz3))#data.frame(orgx)#
Y <- mytot[,,1]

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))

n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 100000; 
burnin <- 25000; 
thinning <- 10

nout <- (iterations-burnin)/thinning
predAPT <- c()

sub <- 42


out_ppmx <- my_dm_ppmx_ct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]), 
                          Xpred = data.frame(X[sub,]), Z = data.frame(Z[-sub,]), 
                          Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                          PPMx = 1, cohesion = 2, alpha = 1, sigma = 0.5,
                          similarity = 2, consim = 2, similparam = vec_par, 
                          calibration = 2, coardegree = 1, modelpriors, update_hierarchy = T,
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

#posterior predictive probabilities ----
A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)



