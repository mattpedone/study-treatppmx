#load all rep1.rda.... and rbind them in case2PPMXUT

#case2PPMXUT<-#utpredAPT_all;
my.pick<-1:K;
wk<-c(0,40,100);
mywk1<-myprob[[1]]%*%wk;  
mywk2<-myprob[[2]]%*%wk;
optrt<-as.numeric( mywk2> mywk1)+1; 
ut.sum<-sum(abs(mywk2-mywk1));ut.diff<- abs(as.numeric(mywk2- mywk1));

#MOT
PPMXCT<-  apply(abs((case2PPMXUT[,3, my.pick]-optrt)),2,sum) 
MOT <- c(round(mean(PPMXCT)), round(sd(PPMXCT), 1))

#MTUg
#MTUg<-cbind(MTU=round(apply(my.result,2,mean)/ut.sum,4),
#SD=round(sqrt(apply(my.result/ut.sum,2,var)),2));
PPMXpp<-  -(2*apply(abs((case2PPMXUT[,3, my.pick]-optrt))*ut.diff,2,sum)-ut.sum);
MTUg <- c(round(mean(PPMXpp/ut.sum), 4), round(sd(PPMXpp/ut.sum), 4))
#NPC
PPMXCUT<-as.vector(countUT(case2PPMXUT));
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))
resPPMX <- rbind(MOT, MTUg,NPC)
colnames(resPPMX) <- c("mean", "sd")
resPPMX
#save(resPPMX,file="resPPMX.rda")

###10
#mean     sd
#MOT  36.0000 3.5000
#MTUg  0.5179 0.0493
#NPC  63.5000 3.5355