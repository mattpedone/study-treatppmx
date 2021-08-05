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
utpred1APT.all<-array(0,dim=c(npat,19,K))

vecadp <- c(1, 2, 10)
idxsc <- 1

for(k in 1:K){
  utpred1APT<-matrix(1,nrow= npat,ncol=19);  ### ut1,ut2,trt,cluster
 for(sub in 1:npat){
   X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
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
   iterations <- 100#0000; 
   burnin <- 0#50000; 
   thinning <- 1#0
   
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
   #utpred1APT[sub,10:12]<-probpre1APT;utpred1APT[sub,13:15]<-probpre2APT;
   #utpred1APT[sub,16:18]<-trtpfySub/sum(trtpfySub);
   #utpred1APT[sub,19]<-max.clus[sub]
   #optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
   #predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
   
   #treatmente prediction with utility function ----
   optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
   cat("optrt: ", optrt, "\n")
   predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
   cat("predtrt: ", predtrt, "\n")
   print(as.numeric(optrt==predtrt))
 }
  utpred1APT.all[,,k]<-utpred1APT
}

case2PPMXUT<-utpred1APT.all;
my.pick<-1:K;
wk<-c(0,40,100);
mywk1<-myprob[[1]]%*%wk;  
#mywk1 <- mywk1[1:50]
mywk2<-myprob[[2]]%*%wk;
#mywk2 <- mywk2[1:50]
optrt<-as.numeric( mywk2> mywk1)+1; 
#optrt <- optrt[1:50]
ut.sum<-sum(abs(mywk2-mywk1));ut.diff<- abs(as.numeric(mywk2- mywk1));
#HCppcont<- apply(abs((case2PPMXUT[,3, my.pick]-optrt)),1,sum);
#MOT
PPMXCT<-  apply(abs((case2PPMXUT[,3, my.pick]-optrt)),2,sum) 
MOT <- c(round(mean(PPMXCT)), round(sd(PPMXCT), 1))

#MTUg
#MTUg<-cbind(MTU=round(apply(my.result,2,mean)/ut.sum,4),
#SD=round(sqrt(apply(my.result/ut.sum,2,var)),2));
PPMXpp<-  -(2*apply(abs((case2PPMXUT[,3, my.pick]-optrt))*ut.diff,2,sum)-ut.sum);
MTUg <- c(round(mean(PPMXpp/ut.sum), 4), round(sd(PPMXpp/ut.sum), 4))
#NPC
PPMXCUT<-as.vector(countUT(case2PPMXUT));
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))
resPPMX <- rbind(MOT, MTUg,NPC)
colnames(resPPMX) <- c("mean", "sd")
resPPMX
#save(resPPMX,file="resPPMX.rda")

###10
#mean     sd
#MOT  36.0000 3.5000
#MTUg  0.5179 0.0493
#NPC  63.5000 3.5355