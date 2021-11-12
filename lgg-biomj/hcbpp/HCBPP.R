## Analyze the data with continuous covariates and gene expressions              #####
## Based on matched data with 158 patients and 79 treated                        #####
## For the implement of HCBPP
######################################################################################

     rm(list=ls()); set.seed(123456);
## R library  
     library("glmnetcr"); 
     library("mvtnorm");
     library(ConsensusClusterPlus);

## include R code for estimation
     source("lgg-biomj/hcbpp/Functions.R");

## load the Example data  

     load("data/LGGdata.rda");
     gene.normAPT<-t(scale(matchRTComp[,16:38]));  ## Predictive features  
             Rapp<-scale(matchRTComp[,c(11,13)]);  ## Prognostic features
###################################################################################
### define outcome, treatment and prognostic variables ################################
         trtAPT<-as.numeric(paste(matchRTComp[,10]));
      outcomAPT<-as.numeric(paste(matchRTComp[,9]));
  
### start the calculation ###########################################################

            wk<-c(0,41,100);                        ## Utility weights
        prior1<-c(1/3,1/3,1/3);                     ## intial prior for non-tagred 
        prior2<-c(1/3,1/3,1/3);                     ## intial prior for tagred 

### hyperparameter for inverse wishart

        kappa0<-1                                   ## kappa for the precision
           mu0<-c(0,0)                              ## assume the mean is know as 0
       n.mysub<-length(trtAPT);

### define some matrix for the calculation
    utpred1APT<-matrix(1,nrow= n.mysub,ncol=19);      ## ut1,ut2,trt,cluster
        HC.sum<-matrix(0,nrow<-n.mysub,ncol=14);

### similarity matrix for all patients
   rst.hcSub<-ConsensusClusterPlus(gene.normAPT,maxK=15,reps=500,pItem=0.90,pFeature=1,
          clusterAlg="hc",distance="pearson",seed=126);
######################################################################################
 for (mysub in 1:n.mysub){
        gene.norm<-gene.normAPT[,-mysub];
              trt<-trtAPT[-mysub];
           outcom<-outcomAPT[-mysub];
     select.sub.n<-length(outcom);
            myRapp<-Rapp[-mysub,];
#################################################################################
       d=gene.norm;
rst.hc<-ConsensusClusterPlus(d,maxK=15,reps=500,pItem=0.90,pFeature=1,
          clusterAlg="hc",distance="pearson",seed=126);
## transform the outcome variables into matrix for easy computation.
   myy<-matrix(0,nrow=select.sub.n,ncol=3);
     for (j in 1:select.sub.n){myy[j,outcom[j]+1]=1}

################################################################################

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
  
########################################################################################
  HC.sum[mysub,]<-c(PreUtBPP(hc2),PreUtBPP(hc3),PreUtBPP(hc4),PreUtBPP(hc5),PreUtBPP(hc6),
PreUtBPP(hc7),PreUtBPP(hc8),PreUtBPP(hc9),PreUtBPP(hc10),PreUtBPP(hc11),PreUtBPP(hc12),
PreUtBPP(hc13),PreUtBPP(hc14),PreUtBPP(hc15));

### determine the rank with the largest summary measures (if same, pick small rank)
           sumtemp<-HC.sum[mysub,1];max.clus<-2;
    for (k in 2:14){
          if (HC.sum[mysub,k]>sumtemp)
             { sumtemp<-HC.sum[mysub,k]; max.clus<-k+1;}
    } 
    
### calculate the predicted probability using the model with largest sum measure
 ## calculate the probablility from the prognostic
      mycovXSub<-Rapp[mysub,];                 ## observed covariates 
        trtcSub<- cbind(outcom,myRapp);
         
         trtc0Sub<-subset(trtcSub,trtcSub[,1]==0); 
      trtc0dstSub<-mymultt(trtc0Sub[,-1],mycovXSub);
         trtc1Sub<-subset(trtcSub,trtcSub[,1]==1); 
      trtc1dstSub<-mymultt(trtc1Sub[,-1],mycovXSub);
         trtc2Sub<-subset(trtcSub,trtcSub[,1]==2); 
      trtc2dstSub<-mymultt(trtc2Sub[,-1],mycovXSub);
        trtpfySub<-c(trtc0dstSub,trtc1dstSub,trtc2dstSub);   ##prob with cov

  myyAPT<-matrix(0,nrow=n.mysub,ncol=3);
     for (m in 1:n.mysub){myyAPT[m,outcomAPT[m]+1]=1};

     mycons1APT<-rst.hcSub[[max.clus]][["consensusMatrix"]];

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
        utpred1APT[mysub,19]<-max.clus
   
   }

#######################################################################################
           mth<-utpred1APT;
           myresults<-cbind(mth,trtAPT+1,as.numeric(outcomAPT>=2))
           pred1<-subset( myresults,myresults[,3]==1);
           table1<-table(pred1[,21],pred1[,20]);
           pred2<-subset( myresults,myresults[,3]==2);
           table2<-table(pred2[,21],pred2[,20]);
             p1<-sum(table1)/(sum(table1)+sum(table2));
             p2<-sum(table2)/(sum(table1)+sum(table2));
########################################################################################
## set prob as 0 for cases that none selected patients had response (all 0s);
## set prob as 1 for cases that none selected patients had non-response (all 1s);
########################################################################################
            if (length(table1)==4){ crt1<-table1[2,1]/sum(table1[,1])};
            if (length(table1)<4) { crt1<-as.numeric(row.names(table1))};

            if (length(table2)==4){  crt2<-table2[2,2]/sum(table2[,2])};
            if (length(table2)<4) {  crt2<-as.numeric(row.names(table2))};
### summary meaures
            return<-crt1*p1+crt2*p2-sum(as.numeric(outcomAPT>=2))/158;
            round(return,3); 

########################################################################################

############################# TRINARY #####################################################
## Calculate CPO counts
   myresults<-cbind(mth,trtAPT+1,as.numeric(outcomAPT));
   ResTable1<-myresults[80:158,c(4:6,21)];              ## received trt1
   ResTable2<-myresults[1:79,c(7:9,21)];                ## received trt2

   ResTable3<-rbind(ResTable1,ResTable2); 

### predicted outcomes;
   MoReult<-as.numeric(apply(ResTable3[,1:3],1,which.max))-1;
   Cresult<-cbind(ResTable3[,4], MoReult);

   table(Cresult[,2],Cresult[,1]);
###  CPO Results 
  sum(diag(table(Cresult[,2],Cresult[,1])))

##########################################################################################
