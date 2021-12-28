## Analysis the data using CONSENSUS with HC                                     #####
## Based on matched data with simulated 152 patients                             #####
## use summary meause in LOOCV for the estimation                                #####
## Junsheng Ma, 5/9/2016                                                         #####
######################################################################################

rm(list=ls()); set.seed(123456);
library("gtools");  library("xtable"); library("mvtnorm");
library("glmnetcr");library(ConsensusClusterPlus);

load("data/scenalt1.RData");     
load("R/ma-scripts/test.RData");  

source("src/countUTtt.R");  

#################################################################################
#### function for covariates## for calculation of the density;
mymultt<-function(Xtrain,X.pred){
  myln<-length(Xtrain[,1]);  myls<-Xtrain;
  mylmu<-apply(myls,2,mean);
  mymun<-myln/(myln+kappa0)*mylmu;
  myS<-cov(myls)*(myln-1); 
  ## for the covariates
  myd<-length(Xtrain[1,])       ## numer of covariates
  nu0<-length(Xtrain[1,])+1       ## numer of covariates +1;
  lambda0<-diag(myd)                ## identity matrix
  #lambda0<-myS*nu0/(myln-1);
  
  kappan<-kappa0+ myln;
  nun<-nu0+  myln;
  lambdn<-lambda0+myS+kappa0*myln/(kappa0+myln)*(mylmu-mu0)%*%t(mylmu-mu0);
  return2<-dmvt(x=X.pred,sigma=(kappan+1)*lambdn/(kappan*(nun-myd+1)),df=nun-myd+1,log=FALSE);
}

################################ setup Parameters ########################################
wk <- c(0,40,100)
prior1 <- c(1/3,1/3,1/3)
prior2<-c(1/3,1/3,1/3)
kappa0<-1
mu0<-c(0,0)
#gene.normAPT<-t(mydata)
#gene.normAPT <- gene.normAPT[,101:152]
#Rapp<-t(rbind(myz2,myz3))
#Rapp <- Rapp[101:152,]
#trtAPT<-as.numeric(trtsgn)-1
#trtAPT <- trtAPT[101:152]
#n<-length(trtAPT);
#trtsgn <- trtsgn[101:152]

d <- scenalt1$pred[[1]]
n <- 28#dim(d)[1]
nrep<-2#30;

utpred1APT.all<-array(0,dim=c(n,19,nrep))

### clustering using CONSENSUS MATRIX method ###################################

rst.hc<-ConsensusClusterPlus(t(d),maxK=15,reps=500,pItem=0.90,pFeature=1,
                             clusterAlg="hc",distance="pearson", 
                             #clusterAlg="km",distance="euclidean", 
                             #clusterAlg="pam",distance="manhattan", 
                             seed=126);

############Analyze the 100 replications ####################################;

for(myrep in 1:nrep)   {  
  #myseed<-123*myrep;    set.seed(myseed);
  
  trtAPT <- scenalt1$trtsgn[[myrep]][125:152]-1
  Rapp <- scenalt1$prog[[myrep]][125:152,]
  #trtAPT <- trtAPT[101:152]
  outcomAPT <- scenalt1$yord[[myrep]][125:152]-1
  #outcomAPT <- outcomAPT[101:152]
  trtsgn <- scenalt1$trtsgn[[myrep]][125:152]
  myoutot <- scenalt1$ymat[[myrep]][125:152,]
  utpred1APT<-matrix(1,nrow= n,ncol=19);  ### ut1,ut2,trt,cluster
  ######################## Measure similarities ###################################
  ### pick the rank with the largest summary measure
  max.clus<-apply(HC.sum.all[,,myrep],1,which.max)+1;
  max.clus <- median(max.clus)
  ##################################################################################### 
  for (mysub in 1:n){
    trt<-trtAPT[-mysub];
    outcom<-outcomAPT[-mysub];
    select.sub.n<-length(outcom);
    myRapp<-Rapp[-mysub,];

    #################################################################################
    myyAPT<-matrix(0,nrow=n,ncol=3);
    for (m in 1:n){myyAPT[m,outcomAPT[m]+1]=1};
    
    mycons1APT<-rst.hc[[max.clus]][["consensusMatrix"]];
    
    ### calculate the predicted probability using the model with largest sum measure
    ## calculate the density/probablility from the prognostic features
    #for(mysub in 1: n){
    mycovXSub<-Rapp[mysub,];                 ## observed covariates 
    trtcSub<- cbind(outcom,myRapp);
    
    trtc0Sub<-subset(trtcSub,trtcSub[,1]==0); 
    trtc0dstSub<-mymultt(trtc0Sub[,-1],mycovXSub);
    trtc1Sub<-subset(trtcSub,trtcSub[,1]==1); 
    trtc1dstSub<-mymultt(trtc1Sub[,-1],mycovXSub);
    trtc2Sub<-subset(trtcSub,trtcSub[,1]==2); 
    trtc2dstSub<-mymultt(trtc2Sub[,-1],mycovXSub);
    trtpfySub<-c(trtc0dstSub,trtc1dstSub,trtc2dstSub);   ##prob with cov
    
    ### calculate the probablility from the predictive features
    myyAPT<-matrix(0,nrow=n,ncol=3);
    for (m in 1:n){myyAPT[m,outcomAPT[m]+1]=1};
    
    #mycons1APT<-rst.hc[[max.clus]][["consensusMatrix"]];
    
    
    totutAPT<-myyAPT*mycons1APT[mysub,][125:152];               ### utility 
    totpreAPT<-totutAPT[-mysub,];  trtpreAPT<-trtAPT[-mysub];
    mytemptAPT<-cbind(totpreAPT,trtpreAPT);
    
    ## calculate alpha hat
    trt1gAPT<-subset(mytemptAPT[,1:3],mytemptAPT[,4]==0);
    trt2gAPT<-subset(mytemptAPT[,1:3],mytemptAPT[,4]==1);
    alphatrt1APT<-apply(trt1gAPT,2,sum)+prior1; 
    probpre1APT<-alphatrt1APT/sum(alphatrt1APT);
    
    alphatrt2APT<-apply(trt2gAPT,2,sum)+prior2; 
    probpre2APT<-alphatrt2APT/sum(alphatrt2APT);
    
    probpre1<- probpre1APT*trtpfySub/sum(probpre1APT*trtpfySub);
    probpre2<- probpre2APT*trtpfySub/sum(probpre2APT*trtpfySub);
    
    ## calculate the utility
    ut1preAPT<-probpre1%*%wk; ut2preAPT<-probpre2%*%wk;
    utpred1APT[mysub,1]<-ut1preAPT;  utpred1APT[mysub,2]<-ut2preAPT;
    if (ut2preAPT>ut1preAPT){utpred1APT[mysub,3]=2};
    utpred1APT[mysub,4:6]<- probpre1;utpred1APT[mysub,7:9]<-probpre2;
    utpred1APT[mysub,10:12]<-probpre1APT;utpred1APT[mysub,13:15]<-probpre2APT;
    utpred1APT[mysub,16:18]<-trtpfySub/sum(trtpfySub);
    utpred1APT[mysub,19]<-max.clus#[mysub]
  }
  
  utpred1APT.all[,,myrep]<-utpred1APT
  
}                        

#####################Save the results#######################################

case2HCppUT<-utpred1APT.all;
my.pick<-1:nrep;
MOT <- MTUg <- NPC <- c()
for(my.pick in 1:nrep){
wk<-c(0,40,100);
myprob <- scenalt1$prob[[my.pick]]

mywk1<-myprob[[1]]%*%wk;  
mywk1 <- mywk1[125:152]
mywk2<-myprob[[2]]%*%wk;
mywk2 <- mywk2[125:152]
optrt<-as.numeric( mywk2> mywk1)+1; 
#optrt <- optrt#[101:152]
ut.sum<-sum(abs(mywk2-mywk1));ut.diff<- abs(as.numeric(mywk2- mywk1));
HCppcont<- sum(abs((case2HCppUT[,3, my.pick]-optrt)));
#MOT
HCppCT<-  sum(abs((case2HCppUT[,3, my.pick]-optrt))) 
MOT[my.pick] <- HCppCT

#MTUg
#MTUg<-cbind(MTU=round(apply(my.result,2,mean)/ut.sum,4),
#SD=round(sqrt(apply(my.result/ut.sum,2,var)),2));
HCpp<-  -(2*sum(abs((case2HCppUT[,3, my.pick]-optrt))*ut.diff)-ut.sum);
MTUg[my.pick] <- HCpp/ut.sum#c(round(mean(HCpp/ut.sum), 4), round(sd(HCpp/ut.sum), 4))
#NPC
outcomAPT <- scenalt1$yord[[my.pick]][125:152]-1
HCppCUT<-as.vector(countUT(case2HCppUT, outcomAPT));
#myoutot troppo lungo
#tutta la roba che importo dovrà essere ridimensionata o comunque divisa in test e validation
NPC[my.pick] <- HCppCUT
}
#resHCpp <- rbind(MOT, MTUg,NPC)

#colnames(resHCpp) <- c("mean", "sd")

#mtug <- HCpp/ut.sum
#save(resHCpp, file="output/simulation-scenarios/scen4/res_ma_hc.rda")
#save(HCppCT, file="output/simulation-scenarios/scen4/mot_ma_hc.rda")
#save(mtug, file="output/simulation-scenarios/scen4/mtug_ma_hc.rda")
#save(HCppCUT, file="output/simulation-scenarios/scen4/npc_ma_hc.rda")


