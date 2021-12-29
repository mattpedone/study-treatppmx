## Analysis the data using CONSENSUS with HC                                    

rm(list=ls())
set.seed(121)
library("gtools")
library("xtable")
library("mvtnorm")
library("glmnetcr")
library(ConsensusClusterPlus)

load("data/scenalt2.RData")    
load("output/simulation-scenarios/train-test/scen-alt-2/ma_hc_tt.RData");  

################################ Functions ########################################
mymultt <- function(Xtrain, X.pred){
  myln <- length(Xtrain[,1])
  myls <- Xtrain
  mylmu <- apply(myls, 2, mean)
  mymun <- myln/(myln+kappa0)*mylmu
  myS <- cov(myls)*(myln-1)
  
  ## for the covariates
  myd <- length(Xtrain[1,])       ## numer of covariates
  nu0 <- length(Xtrain[1,])+1       ## numer of covariates +1;
  lambda0 <- diag(myd)                ## identity matrix
  
  kappan <- kappa0 + myln
  nun <- nu0 + myln
  lambdn <- lambda0 + myS + kappa0*myln/(kappa0 + myln)*(mylmu - mu0)%*%t(mylmu - mu0)
  return2 <- dmvt(x = X.pred, sigma=(kappan + 1)*lambdn/(kappan*(nun - myd + 1)), df = nun - myd + 1, log = FALSE)
}

# calculate the NPC
countUT <- function(resultsum, myoutot){
  myctut <- array(0, dim = c(3, 3, 100))
  myctutSum <- NULL
  for(i in 1:length(my.pick)){
    mycurdata <- resultsum[,,i]
    mypre <- NULL
    pretrt1 <- apply(mycurdata[,4:6], 1, which.max)
    pretrt2 <- apply(mycurdata[,7:9], 1, which.max)
    mypreTall <- cbind(pretrt1, pretrt2)
    for(j in 1:length(trtsgn)){
      mypre[j] <- mypreTall[j, trtsgn[j]]
      }
    sts <- table(mypre, myoutot)
    mysdls <- as.numeric(rownames(sts))
    str1 <- matrix(0, nrow = 3, ncol = 3)
    str1[mysdls,] <- sts
    
    myctut[,,i] <- str1*diag(3)
    myctutSum[i] <- sum(str1*diag(3))
  }
  return <- cbind(myctutSum)
} 

################################ setup Parameters ########################################
wk <- c(0,40,100)
prior1 <- prior2 <- c(1/3,1/3,1/3)
kappa0 <- 1
mu0 <- c(0, 0)
d <- scenalt2$pred[[1]]
n <- 28#dim(d)[1]
nrep <- 30

utpred1APT.all <- array(0, dim = c(n, 19, nrep))

### clustering using CONSENSUS MATRIX method ###################################

rst.hc<-ConsensusClusterPlus(t(d),maxK=15,reps=500,pItem=0.90,pFeature=1,
                             clusterAlg="hc",distance="pearson", 
                             #clusterAlg="km",distance="euclidean", 
                             #clusterAlg="pam",distance="manhattan", 
                             seed=126);


for(myrep in 1:nrep){  
  trtAPT <- scenalt2$trtsgn[[myrep]][125:152]-1
  Rapp <- scenalt2$prog[[myrep]][125:152,]
  outcomAPT <- scenalt2$yord[[myrep]][125:152]-1
  trtsgn <- scenalt2$trtsgn[[myrep]][125:152]
  myoutot <- scenalt2$ymat[[myrep]][125:152,]
  utpred1APT<-matrix(1,nrow= n,ncol=19)  ### ut1,ut2,trt,cluster
  
  ### pick the median rank with the largest summary measure
  max.clus<-apply(HC.sum.all[,,myrep],1,which.max)+1
  max.clus <- median(max.clus)
  
  for (mysub in 1:n){
    trt <- trtAPT[-mysub]
    outcom <- outcomAPT[-mysub]
    select.sub.n <- length(outcom) 
    myRapp <- Rapp[-mysub,]

    mycovXSub <- Rapp[mysub,]                 ## observed covariates 
    trtcSub <- cbind(outcom, myRapp)
    
    trtc0Sub <- subset(trtcSub, trtcSub[,1] == 0)
    trtc0dstSub <- mymultt(trtc0Sub[,-1], mycovXSub)
    trtc1Sub <- subset(trtcSub, trtcSub[,1] == 1)
    trtc1dstSub <- mymultt(trtc1Sub[,-1], mycovXSub)
    trtc2Sub <- subset(trtcSub, trtcSub[,1] == 2)
    trtc2dstSub <- mymultt(trtc2Sub[,-1], mycovXSub)
    trtpfySub <- c(trtc0dstSub, trtc1dstSub, trtc2dstSub)   ##prob with cov
    
    myyAPT <- matrix(0, nrow = n, ncol = 3)
    for(m in 1:n){
      myyAPT[m, outcomAPT[m]+1] = 1
      }
    
    mycons1APT <- rst.hc[[max.clus]][["consensusMatrix"]]
    
    totutAPT <- myyAPT*mycons1APT[124 + mysub, 125:152]               ### utility 
    totpreAPT <- totutAPT[-mysub,]
    trtpreAPT <- trtAPT[-mysub]
    mytemptAPT <- cbind(totpreAPT, trtpreAPT)
    
    ## calculate alpha hat
    trt1gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 0)
    trt2gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 1)
    alphatrt1APT <- apply(trt1gAPT, 2, sum) + prior1
    probpre1APT <- alphatrt1APT/sum(alphatrt1APT)
    
    alphatrt2APT <- apply(trt2gAPT, 2, sum) + prior2 
    probpre2APT <- alphatrt2APT/sum(alphatrt2APT) 
    
    probpre1 <- probpre1APT*trtpfySub/sum(probpre1APT*trtpfySub)
    probpre2 <- probpre2APT*trtpfySub/sum(probpre2APT*trtpfySub)
    
    ## calculate the utility
    ut1preAPT <- probpre1%*%wk
    ut2preAPT <- probpre2%*%wk
    utpred1APT[mysub, 1] <- ut1preAPT
    utpred1APT[mysub, 2] <- ut2preAPT
    if(ut2preAPT > ut1preAPT){
      utpred1APT[mysub, 3] = 2
      }
    utpred1APT[mysub, 4:6] <- probpre1
    utpred1APT[mysub, 7:9] <- probpre2
    utpred1APT[mysub, 10:12] <- probpre1APT
    utpred1APT[mysub, 13:15] <- probpre2APT
    utpred1APT[mysub, 16:18] <- trtpfySub/sum(trtpfySub)
    utpred1APT[mysub, 19] <- max.clus#[mysub]
  }
  
  utpred1APT.all[,,myrep] <- utpred1APT
  
}                        

#####################Save the results#######################################

case2HCppUT <- utpred1APT.all
MOT <- MTUg <- NPC <- c()
for(my.pick in 1:nrep){
  wk<-c(0,40,100)
  myprob <- scenalt2$prob[[my.pick]]
  
  mywk1 <- myprob[[1]]%*%wk
  mywk1 <- mywk1[125:152]
  mywk2 <- myprob[[2]]%*%wk
  mywk2 <- mywk2[125:152]
  optrt <- as.numeric(mywk2 > mywk1) + 1
  
  ut.sum <- sum(abs(mywk2 - mywk1))
  ut.diff <- abs(as.numeric(mywk2 - mywk1))
  HCppcont <- sum(abs((case2HCppUT[,3, my.pick] - optrt)))
  
  MOT[my.pick] <- sum(abs((case2HCppUT[,3, my.pick]-optrt))) 
  MTUg[my.pick] <- (-(2*sum(abs((case2HCppUT[,3, my.pick]-optrt))*ut.diff)-ut.sum))/ut.sum
  outcomAPT <- scenalt2$yord[[my.pick]][125:152] - 1
  HCppCUT <- as.vector(countUT(case2HCppUT, outcomAPT))
  NPC[my.pick] <- HCppCUT
}

save(MOT, file="output/simulation-scenarios/train-test/scen-alt-2/mot_hc.RData")
save(MTUg, file="output/simulation-scenarios/train-test/scen-alt-2/mtug_hc.RData")
save(NPC, file="output/simulation-scenarios/train-test/scen-alt-2/npc_hc.RData")


