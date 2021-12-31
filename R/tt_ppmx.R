rm(list=ls())
set.seed(121)
load("data/scenalt3.RData")
library(treatppmx)
#library(parallel)
#library(doParallel)
library(mcclust)
library(mcclust.ext)

npc2 <- function(output, trtsgn, myoutot){
  K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- array(0, dim = c(3, 3, K))
  myctutSum <- NULL
  for (i in 1:K) {
    mycurdata <- output[, , i]
    mypre <- NULL
    pretrt1 <- apply(mycurdata[, 1:3], 1, which.max)
    pretrt2 <- apply(mycurdata[, 4:6], 1, which.max)
    mypreTall <- cbind(pretrt1, pretrt2)
    for (j in 1:n) {
      mypre[j] <- mypreTall[j, trtsgn[j]]
    }
    sts <- table(mypre, myoutot)
    mysdls <- as.numeric(rownames(sts))
    str1 <- matrix(0, nrow = 3, ncol = 3)
    str1[mysdls, ] <- sts
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}

K <- 30#repliche
npat_pred <- 28

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)
myres0 <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)

for(k in 1:K){
  #predAPT<-matrix(1,nrow= npat,ncol=10);  ### ut1,ut2,trt,cluster
  #cor_all <- parallel::detectCores()-1#cores to be allocated
  #registerDoParallel(cores = cor_all)
  
  #X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
  X_train <- data.frame(scenalt3$pred[[k]][1:124,])
  Z_train <- data.frame(scenalt3$prog[[k]][1:124,])
  Y_train <- data.frame(scenalt3$ymat[[k]][1:124,])
  
  X_test <- data.frame(scenalt3$pred[[k]][125:152,])
  Z_test <- data.frame(scenalt3$prog[[k]][125:152,])
  Y_test <- data.frame(scenalt3$ymat[[k]][125:152,])
  
  trtsgn_train <- scenalt3$trtsgn[[k]][1:124]
  trtsgn_test <- scenalt3$trtsgn[[k]][125:152]
  
  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_L0 <- diag(10, ncol(Y_train))
  modelpriors$hP0_nu0 <- ncol(Y_train) + 2; modelpriors$hP0_V0 <- diag(1.0, ncol(Y_train))
  
  #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
  vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 100000
  burnin <- 50000
  thinning <- 10
  
  nout <- (iterations-burnin)/thinning
  predAPT <- c()
  
  res0 <- ppmxct(y = data.matrix(Y_train), X = data.frame(X_train), 
                              Xpred = data.frame(X_test), Z = data.frame(Z_train), 
                              Zpred = data.frame(Z_test), asstreat = trtsgn_train, #treatment,
                              PPMx = 1, cohesion = 2, alpha = 10, sigma = 0.25,
                              similarity = 2, consim = 2, similparam = vec_par, 
                              calibration = 2, coardegree = 2, modelpriors, 
                              update_hierarchy = T,
                              hsp = T, iter = iterations, burn = burnin, thin = thinning, 
                              mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, nclu_init = 10)
    
    #number of a cluster, mean, binder &varinf ----
    mc <- apply(res0$nclu, 1, mean)
    trt <- trtsgn_train
    num_treat <- table(trt)
    
    cls1 <- t(as.matrix(res0$label[[1]]))[,c(1:num_treat[1])]
    psm1 <- comp.psm(cls1)
    mc_b1 <- minbinder.ext(psm1)
    mc_vi1 <- minVI(psm1)
    
    cls2 <- t(as.matrix(res0$label[[2]]))[,c(1:num_treat[2])]
    psm2 <- comp.psm(cls2)
    mc_b2 <- minbinder.ext(psm2)
    mc_vi2 <- minVI(psm2)
    
    mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
    mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))
    
  ##treatment prediction with utility function ----
  #cat("errori: ", which(rowSums(myres) == 0), "\n")
  myres0[[k]] <- res0
  myres <- apply(res0$pipred, c(1,2,3), median, na.rm=TRUE)
  myclu <- rbind(mc, mc_b, mc_vi)
  myfit <- c(res0$WAIC, res0$lpml)
  #myres <- myres[complete.cases(myres),]
  #sellines <- as.vector(which(rowSums(myres) != 0))
  #sellines <- 1:npat
  A1 <- myres[,, 1]%*%wk 
  A2 <- myres[,, 2]%*%wk
  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(myres[,, 1], myres[,, 2])
  
  nclust_all[k,] <- c(t(myclu))
  gof_all[k,] <- myfit
  #sellines_all[[k]] <- sellines
}

myprob <- scenalt3$prob[[k]]
mywk1 <- myprob[[1]][125:152,]%*%wk
mywk2 <- myprob[[2]][125:152,]%*%wk
optrt <- as.numeric(mywk2 > mywk1) + 1
utsum <- sum(abs(mywk2 - mywk1)) 
utdiff <- abs(as.numeric(mywk2 - mywk1))

#MOT
PPMXCT <- c()
for(k in 1:K){
  #subset <- sellines_all[[k]]
  PPMXCT[k] <-  sum(abs(predAPT_all[, 3, k] - optrt))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for(k in 1:K){
  #subset <- sellines_all[[k]]
  PPMXpp[k] <- -(2*sum(abs((predAPT_all[, 3, k] - optrt)) * utdiff) - utsum);
}

MTUg <- c(round(mean(PPMXpp/utsum), 4), round(sd(PPMXpp/utsum), 4))

#NPC
PPMXCUT <- c()
for(k in 1:K){
  #subset <- sellines_all[[k]]
  temp <- array(0, dim = c(28, 6, 1))
  temp[,,1] <- predAPT_all[, 4:9,k]
  myoutot <- scenalt3$yord[[k]][125:152,]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
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

PPMXpp <- PPMXpp/utsum

save(resPPMX, file="output/simulation-scenarios/train-test/scen-alt-3/res.RData")
save(cluPPMX, file="output/simulation-scenarios/train-test/scen-alt-3/clu.RData")
save(PPMXCT, file="output/simulation-scenarios/train-test/scen-alt-3/mot.RData")
save(PPMXpp, file="output/simulation-scenarios/train-test/scen-alt-3/mtug.RData")
save(PPMXCUT, file="output/simulation-scenarios/train-test/scen-alt-3/npc.RData")
