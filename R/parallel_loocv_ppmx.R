rm(list=ls())
set.seed(121)
load("data/scenario1.rda")
library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

K <- 30 #repliche
npat <- length(trtsgn)

predAPT_all <- array(0, dim = c(npat, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)
myres0 <- sellines_all <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)

for(k in 1:K){
  #predAPT<-matrix(1,nrow= npat,ncol=10);  ### ut1,ut2,trt,cluster
  cor_all <- parallel::detectCores()-1#cores to be allocated
  registerDoParallel(cores = cor_all)
  
  #X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
  X <- data.frame(mydata)
  Z <- data.frame(cbind(myz2, myz3))#data.frame(orgx)#
  Y <- mytot[,,k]
  
  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
  modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(1.0, ncol(Y))
  
  #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
  vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 50000 
  burnin <- 20000
  thinning <- 10
  
  nout <- (iterations-burnin)/thinning
  predAPT <- c()
  
  myres <- foreach(sub = 1:npat, .combine = rbind) %dopar%
    {
    out_ppmx <- tryCatch(expr = ppmxct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]), 
                              Xpred = data.frame(X[sub,]), Z = data.frame(Z[-sub,]), 
                              Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                              PPMx = 1, cohesion = 2, alpha = 10, sigma = 0.25,
                              similarity = 2, consim = 2, similparam = vec_par, 
                              calibration = 2, coardegree = 2, modelpriors, 
                              update_hierarchy = T,
                              hsp = T, iter = iterations, burn = burnin, thin = thinning, 
                              mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, nclu_init = 10), error = function(e){FALSE})
    
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
    #A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)
    A0 <- c(apply(out_ppmx$pipred, c(1,2,3), median, na.rm=TRUE), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)
    ifelse(is.logical(out_ppmx), return(rep(0, 14)), return(A0))
    }
  
  ##treatment prediction with utility function ----
  #cat("errori: ", which(rowSums(myres) == 0), "\n")
  myres0[[k]] <- myres
  myres <- myres[complete.cases(myres),]
  sellines <- as.vector(which(rowSums(myres) != 0))
  #sellines <- 1:npat
  A1 <- myres[sellines, 1:3]%*%wk 
  A2 <- myres[sellines, 4:6]%*%wk
  predAPT_all[sellines, 1, k] <- A1
  predAPT_all[sellines, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[sellines, 3, k] <- myt
  predAPT_all[sellines, 4:9, k] <- myres[sellines, 1:6]
  
  nclust_all[k,] <- apply(myres[sellines, 7:12], 2, mean)
  gof_all[k,] <- apply(myres[sellines, 13:14], 2, mean)
  sellines_all[[k]] <- sellines
}

mywk1 <- myprob[[1]]%*%wk
mywk2 <- myprob[[2]]%*%wk
optrt <- as.numeric(mywk2 > mywk1) + 1
utsum <- sum(abs(mywk2 - mywk1)) 
utdiff <- abs(as.numeric(mywk2 - mywk1))

#MOT
PPMXCT <- c()
for(k in 1:K){
  subset <- sellines_all[[k]]
  PPMXCT[k] <-  sum(abs(predAPT_all[subset, 3, k] - optrt[subset]))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for(k in 1:K){
  subset <- sellines_all[[k]]
  PPMXpp[k] <- -(2*sum(abs((predAPT_all[subset, 3, k] - optrt[subset])) * utdiff[subset]) - utsum);
}

MTUg <- c(round(mean(PPMXpp/utsum), 4), round(sd(PPMXpp/utsum), 4))

#NPC
PPMXCUT <- c()
for(k in 1:K){
  subset <- sellines_all[[k]]
  temp <- array(0, dim = c(length(subset), 6, 1))
  temp[,,1] <- predAPT_all[subset, 4:9,k]
  PPMXCUT[k] <- npc(temp, trtsgn[subset], myoutot[subset,])
}
PPMXCUT <- as.vector(npc(predAPT_all[, 4:9,], trtsgn, myoutot));
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#NCLU
NC <- c(apply(nclust_all[,1:2], 2, mean), apply(nclust_all[,1:2], 2, sd))
BI <- c(apply(nclust_all[,3:4], 2, mean), apply(nclust_all[,3:4], 2, sd))
VI <- c(apply(nclust_all[,5:6], 2, mean), apply(nclust_all[,5:6], 2, sd))

#FIT
WAIC <- c(mean(gof_all[,1]), sd(gof_all[,1]))
lpml <- c(mean(gof_all[,2]), sd(gof_all[,2]))

#results
resPPMX <- rbind(MOT, MTUg, NPC, WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")
resPPMX

cluPPMX <- rbind(NC, BI, VI)
colnames(cluPPMX) <- c("mean trt 1", "mean trt 2", "sd trt 1", "sd trt 2")
cluPPMX <- cluPPMX[, c(1, 3, 2, 4)]
cluPPMX

#save(resPPMX, file="output/simulation-scenarios/scen1/res.RData")
#save(cluPPMX, file="output/simulation-scenarios/scen1/clu.RData")
#save(MOT, file="output/simulation-scenarios/scen1/mot.RData")
#save(MTUg, file="output/simulation-scenarios/scen1/mtug.RData")
#save(NPC, file="output/simulation-scenarios/scen1/npc.RData")