## Analyze the data with continuous covariates and gene expressions              #####
## Based on matched data with 158 patients and 79 treated                        #####
######################################################################################

     rm(list=ls()); set.seed(123456);

## R library  
     library("glmnetcr"); 
     library("mvtnorm");
     library(ConsensusClusterPlus);

## include R code for estimation
     source("Functions.R");

## load the Example data  

     load("LGGdata.rda");  ## 
     gene.normAPT<-t(scale(matchRTComp[,16:38]));  ## Predictive features  
             Rapp<-scale(matchRTComp[,c(11,13)]);  ## Prognostic features           
        gene.norm<-rbind(gene.normAPT,t(Rapp))

### define outcome variables and treatment variables #################################
   
              trt<-as.numeric(paste(matchRTComp[,10]));
           outcom<-as.numeric(paste(matchRTComp[,9]));

### start the calculation ###########################################################

            wk<-c(0,41,100);                        ## Utility weights
        prior1<-c(1/3,1/3,1/3);                     ## intial prior for non-tagred 
        prior2<-c(1/3,1/3,1/3);                     ## intial prior for tagred
  select.sub.n<-length(trt);                        ## number of subjects

######################## Lasso using glmnetcr package###################################
     myx<-gene.norm; myyout<-as.matrix(outcom,ncol=1);
      nsub<-select.sub.n;

### calculate the predictive probabilities
      lasso.result<-LR(1);
      out.response<-as.numeric(outcom>=2);
            SUB.ID<-c(1:select.sub.n)
  
########################################################################################
 lasso.sum<- PreUtLaso( lasso.result);

### ESM (empirical summary measure)
    lasso.sum

######################################################################################## 
        mth=lasso.result;  
   myresults<-cbind(mth,trt+1,out.response,SUB.ID)

######################################################################################## 
## Calculate CPO counts
## number of patients for which the model correctly predicted their observed outcome. 
   myresults<-cbind(mth,trt+1,as.numeric(outcom),SUB.ID);
   ResTable1<-myresults[80:158,c(4:6,11)];              ## received trt1
   ResTable2<-myresults[1:79,c(7:9,11)];                ## received trt2

   ResTable3<-rbind(ResTable1,ResTable2); 

### predicted outcomes;
   MoReult<-as.numeric(apply(ResTable3[,1:3],1,which.max))-1;
   Cresult<-cbind(ResTable3[,4], MoReult);

   table(Cresult[,2],Cresult[,1]);
                                         
###  CPO Results 
  sum(diag(table(Cresult[,2],Cresult[,1])))
##########################################################################################