rm(list=ls())
load("data/SimuOutsce2.rda")
library(treatppmx)
source("src/countUT.R");  

K <- 2 #repliche
npat <- 152
predAPT_all<-array(0,dim=c(npat,10,K))

idxsc <- 1

for(k in 1:K){
  predAPT<-matrix(1,nrow= npat,ncol=10);  ### ut1,ut2,trt,cluster
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
   
   n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
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
   predAPT[sub,10] <- as.double(time_ppmx[3])
   
   #posterior predictive probabilities ----
   A0 <- apply(out_ppmx$ypred, c(1,2,3), mean);#A0
   
   #treatment prediction with utility function ----
   predAPT[sub,1]<-A0[,,1]%*%wk;  predAPT[sub,2]<-A0[,,2]%*%wk;
   if (predAPT[sub,2]>predAPT[sub,1]){predAPT[sub,3]=2};
   predAPT[sub,4:6]<- A0[,,1];predAPT[sub,7:9]<-A0[,,2];
   cat("sample LOOCV: ", sub, "\n")
 }
  predAPT_all[,,k]<-predAPT
}

save(predAPT_all, file="output/scen50/rep1.rda")
