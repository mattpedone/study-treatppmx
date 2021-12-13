### Calculate the Leave-One-Out Cross-validated results
### Date: Feb 12th/2016 
### Author: Junsheng Ma

#############################################################################################
rm(list=ls())
set.seed(101027)
load("data/scenario2.rda")

library(ConsensusClusterPlus); 
library("mvtnorm");

################################ setup Parameters ########################################

wk <- c(0,40,100)
prior1 <- c(1/3,1/3,1/3)
prior2 <- c(1/3,1/3,1/3)
kappa0 <- 1
mu0 <- c(0,0)
gene.normAPT <- t(mydata)
Rapp<-t(rbind(myz2,myz3)) 
trtAPT <- as.numeric(trtsgn)-1
n.mysub <- length(trtAPT)
nrep <- 30

HC.sum.all<-array(0,dim=c(n.mysub,14,nrep))

## define a function that could be used 
###################################################################################

for(foldNumber in 1:nrep){
  myrep<-foldNumber;
  myseed<-123*myrep;       set.seed(myseed);
  outcomAPT<-as.numeric(myoutot[,foldNumber])-1;
  HC.sum<-matrix(0,nrow=n.mysub,ncol=14);
  
  for (mysub in 1:n.mysub){
    gene.norm<-gene.normAPT[,-mysub];
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
    ### clustering using CONSENSUS MATRIX method ###################################
    d=gene.norm;
    rst.hc<-ConsensusClusterPlus(d,maxK=15,reps=500,pItem=0.90,pFeature=1,
                                 #clusterAlg="hc",distance="pearson", 
                                 #clusterAlg="km",distance="euclidean", 
                                 clusterAlg="pam",distance="manhattan",
                                 seed=126);
    
    myy<-matrix(0,nrow=select.sub.n,ncol=3);
    for (j in 1:select.sub.n){myy[j,outcom[j]+1]=1}
    
    con.cluster<-function(cons){ 
      mycons1<-cons;
      utpred1<-matrix(1,nrow=select.sub.n,ncol=12);  ### ut1,ut2,trt
      
      for (myloop in 1:select.sub.n){
        totut<-myy*mycons1[myloop,];               ### utility 
        totpre<-totut[-myloop,];  trtpre<-trt[-myloop];
        
        myRappLp<-myRapp[-myloop,]
        myoutcom<-outcom[-myloop];
        #################################################################################
        mytempt<-cbind(totpre,trtpre,myoutcom,myRappLp);
        
        ## for genes
        trt1g<-subset(mytempt[,1:3],mytempt[,4]==0);
        trt2g<-subset(mytempt[,1:3],mytempt[,4]==1);
        alphatrt1<-apply(trt1g,2,sum)+prior1; probpre1_temp<-alphatrt1/sum(alphatrt1);
        alphatrt2<-apply(trt2g,2,sum)+prior2; probpre2_temp<-alphatrt2/sum(alphatrt2);
        
        #################################################################################
        #################################################################################
        mycovX<-myRapp[myloop,];                  ## observed covariates 
        
        trtc<-subset(mytempt[,-(1:4)]);
        
        trtc0<-subset(trtc,trtc[,1]==0); trtc0dst<-mymultt(trtc0[,-1],mycovX);
        trtc1<-subset(trtc,trtc[,1]==1); trtc1dst<-mymultt(trtc1[,-1],mycovX);
        trtc2<-subset(trtc,trtc[,1]==2); trtc2dst<-mymultt(trtc2[,-1],mycovX);
        trtpfy<-c(trtc0dst,trtc1dst,trtc2dst);   ##prob with cov
        
        #### predicted probability with the genes and covariates;
        probpre1<-probpre1_temp*trtpfy/sum(probpre1_temp*trtpfy);
        probpre2<-probpre2_temp*trtpfy/sum(probpre2_temp*trtpfy);
        
        ## calculate the utility
        ut1pre<-probpre1%*%wk; ut2pre<-probpre2%*%wk;
        utpred1[myloop,1]<-ut1pre;  utpred1[myloop,2]<-ut2pre;
        if (ut2pre>ut1pre){utpred1[myloop,3]=2};
        utpred1[myloop,4:6]<-probpre1;utpred1[myloop,7:9]<-probpre2;
        utpred1[myloop,10:12]<-trtpfy/(sum(trtpfy))
        ut1pre<-ut2pre<-totut<-totpre<-mytempt<-trt1g<-trt2g<-alphatrt1<-alphatrt2<-NULL;
      }
      return<-utpred1;
    }
    
    #################################################################################
    #################################################################################
    
    hc2<-con.cluster(rst.hc[[2]][["consensusMatrix"]]); 
    hc3<-con.cluster(rst.hc[[3]][["consensusMatrix"]]); 
    hc4<-con.cluster(rst.hc[[4]][["consensusMatrix"]]); 
    hc5<-con.cluster(rst.hc[[5]][["consensusMatrix"]]); 
    hc6<-con.cluster(rst.hc[[6]][["consensusMatrix"]]); 
    hc7<-con.cluster(rst.hc[[7]][["consensusMatrix"]]); 
    hc8<-con.cluster(rst.hc[[8]][["consensusMatrix"]]); 
    hc9<-con.cluster(rst.hc[[9]][["consensusMatrix"]]); 
    hc10<-con.cluster(rst.hc[[10]][["consensusMatrix"]]);
    hc11<-con.cluster(rst.hc[[11]][["consensusMatrix"]]);
    hc12<-con.cluster(rst.hc[[12]][["consensusMatrix"]]); 
    hc13<-con.cluster(rst.hc[[13]][["consensusMatrix"]]); 
    hc14<-con.cluster(rst.hc[[14]][["consensusMatrix"]]); 
    hc15<-con.cluster(rst.hc[[15]][["consensusMatrix"]]); 
    
    out.response<-as.numeric(outcom>1);
    SUB.ID<-c(1:select.sub.n)
    ####################################################################################
    
    PreUt<-function(mth){
      ########################################################################################
      ####trt equal 1 nontargeted; 2 targeted#################################################
      ########################################################################################
      mth<-mth;
      myresults<-cbind(mth,trt+1,out.response,SUB.ID)
      pred1<-subset( myresults,myresults[,3]==1);
      table1<-table(pred1[,14],pred1[,13]);
      pred2<-subset( myresults,myresults[,3]==2);
      table2<-table(pred2[,14],pred2[,13]);
      p1<-sum(table1)/(sum(table1)+sum(table2));
      p2<-sum(table2)/(sum(table1)+sum(table2));
      ########################################################################################
      ## set prob as 0 for cases that none selected patients had response (all 0s);
      ## set prob as 1 for cases that none selected patients had non-response (all 1s);
      ########################################################################################
      if (length(table1)==4){ crt1<-table1[2,1]/sum(table1[,1])};
      if (length(table1)<4) { crt1<-0};
      #### if table 1 is empty then set crt1=0;
      if (length(table2)==4){  crt2<-table2[2,2]/sum(table2[,2])};
      if (length(table2)<4) {  crt2<-0};
      
      #### summary meaures
      return<-crt1*p1+crt2*p2-sum(out.response)/length(out.response);
    }
    
    ########################################################################################
    HC.sum[mysub,]<-c(PreUt(hc2),PreUt(hc3),PreUt(hc4),PreUt(hc5),PreUt(hc6),
                      PreUt(hc7),PreUt(hc8),PreUt(hc9),PreUt(hc10),PreUt(hc11),PreUt(hc12),
                      PreUt(hc13),PreUt(hc14),PreUt(hc15));
    
  }
  ########################################################################################
  myresult2<- HC.sum;
  
  # Send a results message back to the master
  
  results <- list(myresult2=myresult2,foldNumber=foldNumber)
  
  HC.sum.all[,,foldNumber] <- HC.sum
}

save(HC.sum.all, file="output/res_ma_pam_scen2.rda")
