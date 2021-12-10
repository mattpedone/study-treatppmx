## Analysis the data using CONSENSUS with HC                                     #####
## Based on matched data with simulated 152 patients                             #####
## use summary meause in LOOCV for the estimation                                #####
## Junsheng Ma, 5/9/2016                                                         #####
######################################################################################

rm(list=ls()); set.seed(123456);
library("gtools");  library("xtable"); library("mvtnorm");
library("glmnetcr");library(ConsensusClusterPlus);

load("data/scenario1.rda");     
load("output/res_ma_hc_scen1.rda");  

source("src/countUT.R");  
################################ setup Parameters ########################################
wk <- c(0,40,100)
prior1 <- c(1/3,1/3,1/3)
prior2<-c(1/3,1/3,1/3)
kappa0<-1
mu0<-c(0,0)
gene.normAPT<-t(mydata)
Rapp<-t(rbind(myz2,myz3))
#Rapp <- Rapp[1:50,]
trtAPT<-as.numeric(trtsgn)-1
#trtAPT <- trtAPT[1:50]
n.mysub<-length(trtAPT);
#trtsgn <- trtsgn[1:50]
nrep<-30;

utpred1APT.all<-array(0,dim=c(n.mysub,19,nrep))

### clustering using CONSENSUS MATRIX method ###################################

rst.hc<-ConsensusClusterPlus(gene.normAPT,maxK=15,reps=500,pItem=0.90,pFeature=1,
                             clusterAlg="hc",distance="pearson",seed=126);

############Analyze the 100 replications ####################################;

for(myrep in 1:nrep)   {  
  myseed<-123*myrep;    set.seed(myseed);
  
  trtAPT<-as.numeric(trtsgn)-1;
  outcomAPT<-as.numeric(myoutot[,myrep])-1;
  utpred1APT<-matrix(1,nrow= n.mysub,ncol=19);  ### ut1,ut2,trt,cluster
  ######################## Measure similarities ###################################
  ### pick the rank with the largest summary measure
  max.clus<-apply(HC.sum.all[,,myrep],1,which.max)+1;
  
  ##################################################################################### 
  for (mysub in 1:n.mysub){
    trt<-trtAPT[-mysub];
    outcom<-outcomAPT[-mysub];
    select.sub.n<-length(outcom);
    myRapp<-Rapp[-mysub,];
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
    #################################################################################
    myyAPT<-matrix(0,nrow=n.mysub,ncol=3);
    for (m in 1:n.mysub){myyAPT[m,outcomAPT[m]+1]=1};
    
    mycons1APT<-rst.hc[[max.clus[mysub]]][["consensusMatrix"]];
    
    ### calculate the predicted probability using the model with largest sum measure
    ## calculate the density/probablility from the prognostic features
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
    myyAPT<-matrix(0,nrow=n.mysub,ncol=3);
    for (m in 1:n.mysub){myyAPT[m,outcomAPT[m]+1]=1};
    
    mycons1APT<-rst.hc[[max.clus[mysub]]][["consensusMatrix"]];
    
    totutAPT<-myyAPT*mycons1APT[mysub,];               ### utility 
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
    utpred1APT[mysub,19]<-max.clus[mysub]
  }
  
  utpred1APT.all[,,myrep]<-utpred1APT
  
}                        

#####################Save the results#######################################

case2HCppUT<-utpred1APT.all;
my.pick<-1:nrep;
wk<-c(0,40,100);
mywk1<-myprob[[1]]%*%wk;  
#mywk1 <- mywk1[1:50]
mywk2<-myprob[[2]]%*%wk;
#mywk2 <- mywk2[1:50]
optrt<-as.numeric( mywk2> mywk1)+1; 
#optrt <- optrt[1:50]
ut.sum<-sum(abs(mywk2-mywk1));ut.diff<- abs(as.numeric(mywk2- mywk1));
HCppcont<- apply(abs((case2HCppUT[,3, my.pick]-optrt)),1,sum);
#MOT
HCppCT<-  apply(abs((case2HCppUT[,3, my.pick]-optrt)),2,sum) 
MOT <- c(round(mean(HCppCT)), round(sd(HCppCT), 1))

#MTUg
#MTUg<-cbind(MTU=round(apply(my.result,2,mean)/ut.sum,4),
#SD=round(sqrt(apply(my.result/ut.sum,2,var)),2));
HCpp<-  -(2*apply(abs((case2HCppUT[,3, my.pick]-optrt))*ut.diff,2,sum)-ut.sum);
MTUg <- c(round(mean(HCpp/ut.sum), 4), round(sd(HCpp/ut.sum), 4))
#NPC
HCppCUT<-as.vector(countUT(case2HCppUT));
NPC <- c(round(mean(HCppCUT), 4), round(sd(HCppCUT), 4))
resHCpp <- rbind(MOT, MTUg,NPC)

colnames(resHCpp) <- c("mean", "sd")

mtug <- HCpp/ut.sum
save(resHCpp, file="output/simulation-scenarios/scen1/res_ma_hc.rda")
save(HCppCT, file="output/simulation-scenarios/scen1/mot_ma_hc.rda")
save(mtug, file="output/simulation-scenarios/scen1/mtug_ma_hc.rda")
save(HCppCUT, file="output/simulation-scenarios/scen1/npc_ma_hc.rda")


