#rm(list=ls())
#
#library(Rcpp)
#library(RcppArmadillo)
rm(list=ls())
load("data/SimuOutsce2.rda")
library(treatppmx)
source("src/countUT.R");  

K <- 2 #repliche
npat <- 152
utpred1APT.all<-array(0,dim=c(npat,9,K))

vecadp <- c(1, 2, 10)
idxsc <- 1

for(k in 1:K){
  utpred1APT<-matrix(1,nrow= npat,ncol=19);  ### ut1,ut2,trt,cluster
 for(sub in 1:npat){
   X <- data.frame(t(mydata))[, -c(51:92)]#data.frame(mydata)#
   Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
   Y <- mytot[,,k]
   idx <- sub
   wk <- c(0, 40, 100)
   
   treattest <- trtsgn[idx]; trt <- trtsgn[-idx]
   Xtest <- data.frame(X[idx,]); X <- data.frame(X[-idx,])
   Ytest <- data.matrix(Y[idx,]); Y <- data.matrix(Y[-idx,])
   Ztest <- data.frame(Z[idx,]); Z <- data.frame(Z[-idx,])
   nobs <- nrow(Y)
   
   optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1; #optrt
   
   modelpriors <- list()
   modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
   modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))
   
   #alpha_DP <- vecadp[alphadp]
   n_aux <- 5
   vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
   #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
   iterations <- 25000; 
   burnin <- 5000; 
   thinning <- 5
   
   nout <- (iterations-burnin)/thinning
   time_ppmx <- system.time(
     out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                               z = Z, zpred = Ztest, asstreat = trt, #treatment,
                               alpha = 1, CC = n_aux, reuse = 1,
                               PPMx = 1, similarity = 2, consim = 2,  gowtot = 1,
                               alphagow = 5, calibration = 2, coardegree = 2,
                               similparam = vec_par, modelpriors, update_hierarchy = T,
                               iter = iterations, burn = burnin, thin = thinning, hsp = T))
   #risultati[idxsc,9,k] <- as.double(time_ppmx[3])
   
   #posterior predictive probabilities ----
   A0 <- apply(out_ppmx$ypred, c(1,2,3), mean);#A0
   
   #treatmente prediction with utility function ----
   ut1preAPT<-A0[,,1]%*%wk; ut2preAPT<-A0[,,2]%*%wk;
   utpred1APT[sub,1]<-ut1preAPT;  utpred1APT[sub,2]<-ut2preAPT;
   if (ut2preAPT>ut1preAPT){utpred1APT[sub,3]=2};
   utpred1APT[sub,4:6]<- A0[,,1];utpred1APT[sub,7:9]<-A0[,,2];
   
   #treatmente prediction with utility function ----
   optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
   #cat("optrt: ", optrt, "\n")
   predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
   #cat("predtrt: ", predtrt, "\n")
   cat("repl: ", k, "sub: ", sub, "pred: ", as.numeric(optrt==predtrt), "\n")
 }
  utpred1APT.all[,,k]<-utpred1APT
  cat("k: ", k, "\n")
}

save(utpred1APT.all, file="rep1.rda")
