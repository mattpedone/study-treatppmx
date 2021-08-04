### Calculate the Leave-One-Out Cross-validated results
### Date: Feb 12th/2016 
### Author: Junsheng Ma

#############################################################################################
# Initialize MPI
   rm(list=ls());set.seed(101027); getwd();
# Load the R MPI package if it is not already loaded.
if (!is.loaded("mpi_initialize")) {
    library("Rmpi")
    }
###########################################################################
  load("/scratch/biostatistics/jma4/project5JUNE1/SimuCase2/lasso/SimuOutsce2.rda");       

   library("glmnetcr");
################################ setup Parameters ########################################

              wk<-c(0,40,100);
          prior1<-c(1/3,1/3,1/3);   
          prior2<-c(1/3,1/3,1/3); 
     gene.normAPT<-rbind(mydata[1:90,],myx2,myx3);
          trtAPT<-as.numeric(trtsgn)-1;
         n.mysub<-length(trtAPT);
             nrep<-100;
 
################ Start Rmpi part here #############################################################
# Spawn as many slaves as possible
mpi.spawn.Rslaves(nslaves=48)
#mpi.spawn.Rslaves()

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves()
        }
        print("Please use mpi.quit() to quit R")
        .Call("mpi_finalize")
    }
}

# Function the slaves will call to perform a validation on the
# fold equal to their slave number.
# Assumes: thedata,fold,foldNumber,p
foldslave <- function() {
    # Note the use of the tag for sent messages:
   #     1=ready_for_task, 2=done_task, 3=exiting
    # Note the use of the tag for received messages:
    #     1=task, 2=done_tasks
     junk <- 0;     done <- 0
     while (done != 1) {
        # Signal being ready to receive a new task
         mpi.send.Robj(junk,0,1)
        # Receive a task
         
        task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]

        if (tag == 1) {
            foldNumber <- task$foldNumber
            set.seed(foldNumber)

## define a function that could be used 
###################################################################################

           library("glmnetcr");
            myrep<-foldNumber;
           myseed<-123*myrep;       set.seed(myseed);
        outcomAPT<-as.numeric(myoutot[,foldNumber])-1;
              myx<-gene.normAPT;
           myyout<-as.matrix(outcomAPT,ncol=1);
             nsub<-n.mysub;

LR<-function(alpha)  {
     myalpha<-alpha;
     utpredLas<-matrix(1,nrow=nsub,ncol=9);  ### ut1,ut2,trt

    for (myloop in 1:nsub){
         trtpre<-trtAPT[-myloop]+1; ### 2 for targeted treatment
        
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

  ridge.result<-LR(0);  lasso.result<-LR(1);

 ########################################################################################       
  
            # Send a results message back to the master
            
    results <- list(ridge.result=ridge.result,lasso.result=lasso.result,foldNumber=foldNumber)
            mpi.send.Robj(results,0,2)
            }
        else if (tag == 2) {
            done <- 1
            }
        # We'll just ignore any unknown messages
        }
    
mpi.send.Robj(junk,0,3)
    }

# Now, send the data to the slaves
    mpi.bcast.Robj2slave(wk);
    mpi.bcast.Robj2slave(prior1);
    mpi.bcast.Robj2slave(prior2);          
    mpi.bcast.Robj2slave(gene.normAPT);
    mpi.bcast.Robj2slave(trtAPT);
    mpi.bcast.Robj2slave(n.mysub);
    mpi.bcast.Robj2slave(myoutot); 
    mpi.bcast.Robj2slave(nrep);


# send the function to the slaves
mpi.bcast.Robj2slave(foldslave)

# Call the function in all the slaves to get them ready to
# undertake tasks
mpi.bcast.cmd(foldslave())

# Create task list
tasks <- vector('list')
for (i in 1:nrep) {
    tasks[[i]] <- list(foldNumber=i)
    }
# Create data structure to store the results #####################################

    Ridge.all<-Lasso.all<-array(0,dim=c(n.mysub,9,nrep));           ### p-value and RMST
junk <- 0 ;  closed_slaves <- 0 ; n_slaves <- mpi.comm.size()-1
 
while (closed_slaves < n_slaves) {
     # Receive a message from a slave
     message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
     message_info <- mpi.get.sourcetag()
     slave_id <- message_info[1]
     tag <- message_info[2]
    if (tag == 1) {
         # slave is ready for a task. Give it the next task, or tell it tasks
         # are done if there are none.
         if (length(tasks) > 0) {
             # Send a task, and then remove it from the task list
            mpi.send.Robj(tasks[[1]], slave_id, 1);
             tasks[[1]] <- NULL
             }
         else {
             mpi.send.Robj(junk, slave_id, 2)
             }
         }
     else if (tag == 2) {
         # The message contains results. Do something with the results.
         # Store them in the data structure
           mysample<-NULL;

        foldNumber <- message$foldNumber;
          adrs<- paste("/scratch/biostatistics/jma4/project5JUNE/SimuCase2/lasso/",
"Rgresult",foldNumber,".txt",sep="");
           Ridge.all[,,foldNumber]<-message$ridge.result;
           Lasso.all[,,foldNumber]<-message$lasso.result;
          #write.table( Ridge.all[,,foldNumber],file=adrs);

         } 
    else if (tag == 3) {
         # A slave has closed down.
        closed_slaves <- closed_slaves + 1
         }
     } 
## end of Rmpi and calculate final results ##############################################

  save(Ridge.all,Lasso.all,
file="/scratch/biostatistics/jma4/project5JUNE1/SimuCase2/lasso/SimuCase2Lasso.rda")

# Tell all slaves to close down, and exit the program
mpi.close.Rslaves()
mpi.quit(save="no")

