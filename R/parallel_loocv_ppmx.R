rm(list=ls())
load("data/modscenario2.rda")
library(treatppmx)
library(parallel)
library(doParallel)

K <- 3 #repliche
npat <- length(trtsgn)
predAPT_all <- array(0, dim = c(npat, 9, K))

wk <- c(0, 40, 100)

for(k in 1:K){
  #predAPT<-matrix(1,nrow= npat,ncol=10);  ### ut1,ut2,trt,cluster
  cor_all <- parallel::detectCores()-1#cores to be allocated
  registerDoParallel(cores = cor_all)
  
  #X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
  X <- data.frame(mydata)
  Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
  Y <- mytot[,,k]
  
  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
  modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))
  
  n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
  vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 10#0000; 
  burnin <- 1#5000; 
  thinning <- 1#0
  
  nout <- (iterations-burnin)/thinning
  predAPT <- c()
  
  myres <- foreach(sub = 1:npat, .combine = rbind) %dopar%
    {
    out_ppmx <- my_dm_ppmx_ct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]), Xpred = data.frame(X[sub,]),
                                z = data.frame(Z[-sub,]), zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                                alpha = 1, CC = n_aux, reuse = 1,
                                PPMx = 1, similarity = 2, consim = 2,  #gowtot = 1,
                                #alphagow = 5, 
                              calibration = 2, coardegree = 1,
                                similparam = vec_par, modelpriors, update_hierarchy = T,
                                iter = iterations, burn = burnin, thin = thinning, hsp = T)
    #posterior predictive probabilities ----
    A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean));#A0
    return(A0)
    }
  
  ##treatment prediction with utility function ----
  A1 <- myres[,1:3]%*%wk; A2 <- myres[,4:6]%*%wk
  predAPT_all[,1,k]<-A1
  predAPT_all[,2,k]<-A2
  myt <- as.numeric(A1<A2) +1
  predAPT_all[,3,k]<-myt
  predAPT_all[,4:9,k]<-myres
}

#save(predAPT_all, file="output/scen50/rep1.rda")
