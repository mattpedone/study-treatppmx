#### function for covariates## for calculation of the density;
  mymultt<-function(Xtrain,X.pred){
             myln<-length(Xtrain[,1]);  myls<-Xtrain;
             mylmu<-apply(myls,2,mean);
             mymun<-myln/(myln+kappa0)*mylmu;
             myS<-cov(myls)*(myln-1); 
## for the covariates
                myd<-length(Xtrain[1,])       ## numer of covariates
                nu0<-length(Xtrain[1,])+1;    ## numer of covariates +1;
            lambda0<-diag(myd)                ## identity matrix
          #lambda0<-myS*nu0/(myln-1);
        
            kappan<-kappa0+ myln;
               nun<-nu0+  myln;
            lambdn<-lambda0+myS+kappa0*myln/(kappa0+myln)*(mylmu-mu0)%*%t(mylmu-mu0);
return2<-dmvt(x=X.pred,sigma=(kappan+1)*lambdn/(kappan*(nun-myd+1)),df=nun-myd+1,log=FALSE);
   }

#################################################################################
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

####################################################################################
### for BPP 
   PreUtBPP<-function(mth){
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

####################################################################################
### for Lasso 
   PreUtLaso<-function(mth){
########################################################################################
####trt equal 1 nontargeted; 2 targeted#################################################
########################################################################################
           myresults<-cbind(mth,trt+1,out.response,SUB.ID)
           pred1<-subset( myresults,myresults[,3]==1);
           table1<-table(pred1[,11],pred1[,10]);
           pred2<-subset( myresults,myresults[,3]==2);
           table2<-table(pred2[,11],pred2[,10]);
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

#### summary meaures
            return<-crt1*p1+crt2*p2-sum(out.response)/length(out.response);
     }

### function for Ridge and Lasso method, alpha=0 is Ridge, and alpha=1 is Lasso
 LR<-function(alpha)  {
     myalpha<-alpha;
     utpredLas<-matrix(1,nrow=nsub,ncol=9);  ### ut1,ut2,trt

    for (myloop in 1:nsub){
         trtpre<-trt[-myloop]+1; ### 2 for targeted treatment
        
         myxp<-cbind(trtpre,t(myx)[-myloop,]); 
               myxp1<-subset(myxp[,-1],myxp[,1]==1); 
               myxp2<-subset(myxp[,-1],myxp[,1]==2);
         myyoutp<-cbind(trtpre,myyout[-myloop,]);
               myyoutp1<-subset(myyoutp[,-1],myyoutp[,1]==1);
               myyoutp2<-subset(myyoutp[,-1],myyoutp[,1]==2);
 
 ## need to be factors for the outcome variables
         y1<-as.factor(myyoutp1);y2<-as.factor(myyoutp2);

 ################################################################# 
  ##in case we observe only two categories use if to control if
         ncat1<-length(table(y1)); ncat2<-length(table(y2));
         if(ncat1<3 |ncat2<3 ) break; 
 #################################################################
         NewX<-t(myx)[myloop,];
      myfit1<-glmnet.cr(myxp1, y1,method="backward",maxit=2000,alpha=myalpha);
            My.select1<-select.glmnet.cr(myfit1,which="AIC");
           probpre1<-predict(myfit1,newx<-NewX)$probs[,,My.select1];

 #print( nonzero.glmnet.cr(myfit1,s=My.select1))
       myfit2<-glmnet.cr(myxp2, y2,method="backward",maxit=2000,alpha=myalpha);
          My.select2<-select.glmnet.cr(myfit2,which="AIC");
          probpre2<-predict(myfit2,newx<-NewX)$probs[,,My.select2];
 ## print(nonzero.glmnet.cr(myfit2,s=My.select2));
 
 ## calculate the utility
         ut1pre<-probpre1%*%wk; ut2pre<-probpre2%*%wk;
         utpredLas[myloop,1]<-ut1pre;  utpredLas[myloop,2]<-ut2pre;
         if (ut2pre>ut1pre){utpredLas[myloop,3]=2};
          utpredLas[myloop,4:6]<-probpre1;
          utpredLas[myloop,7:9]<-probpre2;
      }

    return<-utpredLas;
  }









