### Calculate the Training-Test results
rm(list=ls())
set.seed(121)

load("data/scenalt4.RData")

library(ConsensusClusterPlus); 
library("mvtnorm");

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

con.cluster <- function(cons, ymat = s_train_ymat,  yvec = s_train_yord, 
                      prog = s_train_prog, trt = s_train_trt){ 
  ni <- dim(cons)[1]
  utpred1 <- matrix(1, nrow = ni, ncol=12)  ## ut1,ut2,trt
  
  for (i in 1:ni){
    totut <- ymat*cons[i,]  ### utility 
    totpre <- totut[-i,]
    trtpre <- trt[-i]
    
    l_train_prog <- prog[-i,]
    l_yvec <- yvec[-i];
    
    mytempt <- cbind(totpre, trtpre, l_yvec, l_train_prog)
    
    ## for genes
    trt1g <- subset(mytempt[,1:3], mytempt[,4] == 0)
    trt2g <- subset(mytempt[,1:3], mytempt[,4] == 1)
    alphatrt1 <- apply(trt1g, 2, sum) + prior1
    probpre1_temp <- alphatrt1/sum(alphatrt1)
    alphatrt2 <- apply(trt2g, 2, sum) + prior2
    probpre2_temp <- alphatrt2/sum(alphatrt2)
    
    mycovX <- prog[i,];                  ## observed covariates 
    
    trtc <- subset(mytempt[,-(1:4)])
    
    trtc0 <- subset(trtc,trtc[,1] == 0)
    trtc0dst <- mymultt(trtc0[,-1], mycovX)
    trtc1 <- subset(trtc,trtc[,1] == 1) 
    trtc1dst <- mymultt(trtc1[,-1], mycovX)
    trtc2 <- subset(trtc, trtc[,1] == 2)
    trtc2dst <- mymultt(trtc2[,-1], mycovX)
    trtpfy <- c(trtc0dst, trtc1dst, trtc2dst)   ##prob with cov
    
    probpre1 <- probpre1_temp*trtpfy/sum(probpre1_temp*trtpfy)
    probpre2 <- probpre2_temp*trtpfy/sum(probpre2_temp*trtpfy)
    
    ## calculate the utility
    ut1pre <- probpre1%*%wk
    ut2pre <- probpre2%*%wk
    utpred1[i, 1] <- ut1pre
    utpred1[i, 2] <- ut2pre
    if(ut2pre > ut1pre){
      utpred1[i, 3] = 2
      }
    utpred1[i, 4:6] <- probpre1
    utpred1[i, 7:9] <- probpre2
    utpred1[i, 10:12] <- trtpfy/(sum(trtpfy))
    ut1pre<-ut2pre<-totut<-totpre<-mytempt<-trt1g<-trt2g<-alphatrt1<-alphatrt2<-NULL;
  }
  return <- utpred1;
}

PreUt<-function(mth){
  #mth<-mth;
  myresults <- cbind(mth, s_train_trt+1, out.response, SUB.ID)
  pred1 <- subset(myresults, myresults[,3] == 1)
  table1 <- table(pred1[,14], pred1[,13])
  pred2 <- subset(myresults, myresults[,3] == 2)
  table2 <- table(pred2[,14], pred2[,13])
  p1 <- sum(table1)/(sum(table1) + sum(table2))
  p2 <- sum(table2)/(sum(table1) + sum(table2))
  ## set prob as 0 for cases that none selected patients had response (all 0s);
  ## set prob as 1 for cases that none selected patients had non-response (all 1s);
  if(length(table1) == 4){
    crt1<-table1[2,1]/sum(table1[,1])
    }
  if(length(table1) < 4){
    crt1 <- 0
    }
  #### if table 1 is empty then set crt1=0;
  if(length(table2) == 4){
    crt2 <- table2[2,2]/sum(table2[,2])
    }
  if(length(table2) < 4){
    crt2 <- 0
    }
  
  #### summary meaures
  return <- crt1*p1 + crt2*p2 - sum(out.response)/length(out.response)
}

################################ setup Parameters ########################################

wk <- c(0,40,100)
prior1 <- c(1/3,1/3,1/3)
prior2 <- c(1/3,1/3,1/3)
kappa0 <- 1
mu0 <- c(0,0)
n <- 124
nrep <- 30

HC.sum.all <- array(0, dim = c(n, 14, nrep))


for(rep in 1:nrep){
  train_pred <- scenalt4$pred[[rep]][1:124,]
  train_prog <- scenalt4$prog[[rep]][1:124,]
  train_yord <- scenalt4$yord[[rep]][1:124]-1
  train_ymat <- scenalt4$ymat[[rep]][1:124,]
  train_trt <- scenalt4$trtsgn[[rep]][1:124]-1
  
  HC.sum<-matrix(0,nrow=n,ncol=14)
  
  for (mysub in 1:n){
    s_train_pred <- train_pred[-mysub,]
    s_train_prog <- train_prog[-mysub,]
    s_train_yord <- train_yord[-mysub]
    s_train_ymat <- train_ymat[-mysub,]
    s_train_trt <- train_trt[-mysub]
    
    ### clustering using CONSENSUS MATRIX method ###################################
    con_clu <- ConsensusClusterPlus(t(s_train_pred),maxK=15,reps=500,pItem=0.90,pFeature=1,
                                    #clusterAlg="hc",distance="pearson", 
                                    clusterAlg="km",distance="euclidean", 
                                    #clusterAlg="pam",distance="manhattan",
                                    seed=126)
    
    hc2 <- con.cluster(con_clu[[2]][["consensusMatrix"]])
    hc3 <- con.cluster(con_clu[[3]][["consensusMatrix"]])
    hc4 <- con.cluster(con_clu[[4]][["consensusMatrix"]])
    hc5 <- con.cluster(con_clu[[5]][["consensusMatrix"]])
    hc6 <- con.cluster(con_clu[[6]][["consensusMatrix"]])
    hc7 <- con.cluster(con_clu[[7]][["consensusMatrix"]])
    hc8 <- con.cluster(con_clu[[8]][["consensusMatrix"]])
    hc9 <- con.cluster(con_clu[[9]][["consensusMatrix"]])
    hc10 <- con.cluster(con_clu[[10]][["consensusMatrix"]])
    hc11 <- con.cluster(con_clu[[11]][["consensusMatrix"]])
    hc12 <- con.cluster(con_clu[[12]][["consensusMatrix"]])
    hc13 <- con.cluster(con_clu[[13]][["consensusMatrix"]])
    hc14 <- con.cluster(con_clu[[14]][["consensusMatrix"]])
    hc15 <- con.cluster(con_clu[[15]][["consensusMatrix"]])
    
    out.response <- as.numeric(s_train_yord > 1)
    SUB.ID <- c(1:(n-1))
    
    HC.sum[mysub,] <- c(PreUt(hc2), PreUt(hc3), PreUt(hc4), PreUt(hc5), PreUt(hc6),
                      PreUt(hc7), PreUt(hc8), PreUt(hc9), PreUt(hc10), PreUt(hc11), 
                      PreUt(hc12), PreUt(hc13), PreUt(hc14), PreUt(hc15))
    
  }
  
  myresult2<- HC.sum
  
  results <- list(myresult2 = myresult2, rep = rep)
  
  HC.sum.all[,,rep] <- HC.sum
}

save(HC.sum.all, file = "output/simulation-scenarios/train-test/scen-alt-4/ma_km_tt.RData")

