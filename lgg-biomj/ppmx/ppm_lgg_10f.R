rm(list=ls())
set.seed(121)

library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)
library(doRNG)

load("data/LGGdata.rda")
#name <- c("v01")

matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
#trtsgn <- ifelse(trtsgn == 1, 2, 1)
npat <- length(trtsgn)
#trtsgn <- sample(1:2, npat, replace = TRUE)

K <- 10#numero di fold

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

registerDoParallel(cores = K)#alloco solo core necessari

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)
#load("/home/matt/Dropbox/PHD/study-treatppmx/output/lgg12aprs121.RData")
#nout <- 1600
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

myres0 <- foreach(k = 1:10) %dopar%
  {
    currfold <- (vectf[k]:(vectf[k+1]-1))
    X_train <- data.frame(X[-currfold,])
    Z_train <- data.frame(Z[-currfold,])
    Y_train <- data.frame(Y[-currfold,])
    
    X_test <- data.frame(X[currfold,])
    Z_test <- data.frame(Z[currfold,])
    Y_test <- data.frame(Y[currfold,])
    
    trtsgn_train <- trtsgn[-currfold]
    trtsgn_test <- trtsgn[currfold]
    
    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_L0 <- diag(10, ncol(Y_train))
    modelpriors$hP0_nu0 <- ncol(Y_train) + 2; modelpriors$hP0_V0 <- diag(1, ncol(Y_train))
    
    #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 4000
    thinning <- 5
    
    nout <- (iterations-burnin)/thinning
    predAPT <- c()
    
    res0 <- tryCatch(expr = ppmxct(y = data.matrix(Y_train), X = data.frame(X_train), 
                                   Xpred = data.frame(X_test), Z = data.frame(Z_train), 
                                   Zpred = data.frame(Z_test), asstreat = trtsgn_train, #treatment,
                                   PPMx = 0, cohesion = 1, kappa = c(.1, 10, 5, 1), sigma = c(0.01, .5, 6),
                                   similarity = 2, consim = 2, similparam = vec_par, 
                                   calibration = 2, coardegree = 2, modelpriors, 
                                   update_hierarchy = T,
                                   hsp = F, iter = iterations, burn = burnin, thin = thinning, 
                                   mhtunepar = c(0.05, 0.05), CC = 3, reuse = 1, 
                                   nclu_init = 10), error = function(e){FALSE})
    return(res0)
  }

for(k in 1:K){
  currfold <- (vectf[k]:(vectf[k+1]-1))
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)
  
  cls1 <- t(as.matrix(res0$label[[1]]))[,c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  mc_b1 <- minbinder.ext(psm1)
  mc_vi1 <- minVI(psm1)
  
  cls2 <- t(as.matrix(res0$label[[2]]))[,c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  mc_b2 <- minbinder.ext(psm2)
  mc_vi2 <- minVI(psm2)
  
  mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
  mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))
  
  myres <- apply(res0$pipred, c(1,2,3), median, na.rm=TRUE)
  myclu <- rbind(mc, mc_b, mc_vi)
  myfit <- c(res0$WAIC, res0$lpml)
  A1 <- myres[,, 1]%*%wk 
  A2 <- myres[,, 2]%*%wk
  predAPT_all[currfold, 1] <- A1
  predAPT_all[currfold, 2] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[currfold, 3] <- myt
  predAPT_all[currfold, 4:9] <- cbind(myres[,, 1], myres[,, 2])
  
  nclust_all[k,] <- c(t(myclu))
  gof_all[k,] <- myfit
  
  #myprob <- simdata$prob[[k]]
}

#NPC
npc_tf <- function(output, trtsgn, myoutot){
  #K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- matrix(0, nrow = 3, ncol = 3)
  myctutSum <- NULL
  mycurdata <- output
  mypre <- NULL
  pretrt1 <- apply(mycurdata[, 1:3], 1, which.max)
  pretrt2 <- apply(mycurdata[, 4:6], 1, which.max)
  mypreTall <- cbind(pretrt1, pretrt2)
  for (j in 1:n) {
    mypre[j] <- mypreTall[j, trtsgn[j]]
    }
  sts <- table(mypre, myoutot)
  mysdls <- as.numeric(rownames(sts))
  str1 <- matrix(0, nrow = 3, ncol = 3)
  str1[mysdls, ] <- sts
  myctut <- str1 * diag(3)
  myctutSum <- sum(str1 * diag(3))
  #res <- cbind(myctutSum)
  return(myctutSum)
}


PPMXCUT <- c()
temp <- matrix(0, nrow = npat, ncol = 6)
for(k in 1:K){
  currfold <- (vectf[k]:(vectf[k+1]-1))
  myoutot <- as.numeric(matchRTComp[currfold,9])#simdata$yord[[k]][131:158,]
  trtsgn_test <- trtsgn[currfold]#simdata$trtsgn[[k]][131:158]
  temp[currfold,] <- predAPT_all[currfold, 4:9]
}
NPC <- npc_tf(temp, trtsgn, as.numeric(matchRTComp[,9]))
#NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#ESM
# non ho definito come respondent anche i partial responent
myoutot <- as.numeric(matchRTComp[,9])#simdata$yord[[k]][131:158,]
mytab <- cbind(myass = predAPT_all[,3], rndass = trtsgn, resp = as.numeric(myoutot>2))
pred1 <- subset(mytab, mytab[,1]==1)
table1 <- table(pred1[,3],pred1[,2])
pred2 <- subset(mytab, mytab[,1]==2)
table2 <- table(pred2[,3], pred2[,2])
p1 <- sum(table1)/(sum(table1)+sum(table2))
p2 <- sum(table2)/(sum(table1)+sum(table2))

if(length(table1) == 4){
  crt1 <- table1[2,1]/sum(table1[,1])
}
if(length(table1) < 4){
  crt1 <- as.numeric(row.names(table1))
}

if(length(table2) == 4){
  crt2 <- table2[2,2]/sum(table2[,2])
}
if(length(table2) < 4){
  crt2 <- as.numeric(row.names(table2))
}
ESM <- c(crt1*p1 + crt2*p2 - sum(as.numeric(myoutot>2))/npat)
### summary meaures

#FIT
WAIC <- c(mean(gof_all[,1]), sd(gof_all[,1]))
lpml <- c(mean(gof_all[,2]), sd(gof_all[,2]))

#results
resPPMX <- rbind(WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")

cluPPMX <- nclust_all[,-c(3,4)]
clu <- apply(cluPPMX, 2, mean)
clu <- rbind(clu, apply(cluPPMX, 2, sd))
colnames(clu) <- c("avg # trt 1", "avg # trt 2", "VI trt 1", "VI trt 2")
NPC; ESM; resPPMX; clu

#save(myres0, file = "output/lgg_analysis_24mar.RData")
nout <- 1600
k <- 7
out_ppmx <- myres0[[k]]

currfold <- (vectf[k]:(vectf[k+1]-1))
#res0 <- myres0[[k]]
#number of a cluster, mean, binder &varinf ----
mc <- apply(out_ppmx$nclu, 1, mean)
trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
num_treat <- table(trt)

# Posterior clustering ----
num_treat <- table(trt)
cls <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
psm <- comp.psm(cls)
mc_b <- minbinder.ext(psm); max(mc_b$cl)
mc_vi <- minVI(psm); max(mc_vi$cl)
reord <- c()
reordb <- c()
for(i in 1:max(mc_vi$cl)){
  reord <- c(reord, which(mc_vi$cl == i))
}
for(i in 1:max(mc_b$cl)){
  reordb <- c(reordb, which(mc_b$cl == i))
}

cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
psm2 <- comp.psm(cls2)
mc_b2 <- minbinder.ext(psm2); max(mc_b2$cl)
mc_vi2 <- minVI(psm2); max(mc_vi2$cl)
reord2 <- c()
reord2b <- c()
for(i in 1:max(mc_vi2$cl)){
  reord2 <- c(reord2, which(mc_vi2$cl == i))
}
for(i in 1:max(mc_b2$cl)){
  reord2b <- c(reord2b, which(mc_b2$cl == i))
}

# Co-occurence plot ----
data <- t(out_ppmx$label[[1]])
data <- data[,reordb]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(coincidences)
c1 <- ggplot(mC, aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  theme_classic() +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 1")

cp1 <- c1 + labs(fill = "Probability")

data <- t(out_ppmx$label[[2]])
data <- data[,reord2b]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
c2 <- ggplot(melt(coincidences), aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  theme_classic() +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 2")

cp2 <- c2 + labs(fill = "Probability")

corrp_vi <- ggpubr::ggarrange(cp1, cp2, nrow=1, ncol = 2, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
corrp_vi
#ggsave(corrp_vi, device = "pdf", path = "figs", filename = "corr_plot.pdf")

# Traceplot for the number of clusters ----
df <- data.frame(t(out_ppmx$nclu))
colnames(df) <- c("treatment 1", "treatment 2")
df <- cbind(Index = as.numeric(row.names(df)), df)
df <- reshape2::melt(df, id.vars="Index")
colnames(df) <- c("Index", "Treatment", "value")
clu_tp <- ggplot2::ggplot(df, aes(x = Index, y = value, col = Treatment)) + 
  geom_line() + theme_classic() + 
  xlab("Iterations") + ylab("# of clusters") 

clu_tp
#ggsave(clu_tp, device = "pdf", path = "figs", filename = "clu_trplot.pdf")

# Posterior frequency for (\kappa, \sigma) ----

df_sigma <- data.frame(table(out_ppmx$sigmangg[1,]))
colnames(df_sigma) <- c("sigma", "frequency")
ms1 <- ggplot(df_sigma, aes(x=sigma, y=frequency)) + 
  geom_segment(aes(x=sigma, xend=sigma, y=0, yend=frequency/nout)) + 
  ylab("proportion") + xlab(expression(sigma))
df_sigma <- data.frame(table(out_ppmx$sigmangg[2,]))
colnames(df_sigma) <- c("sigma", "frequency")
ms2 <- ggplot(df_sigma, aes(x=sigma, y=frequency)) + 
  geom_segment(aes(x=sigma, xend=sigma, y=0, yend=frequency/nout))+ 
  ylab("proportion")+ xlab(expression(sigma))

df_kappa <- data.frame(table(out_ppmx$kappangg[1,]))
colnames(df_kappa) <- c("kappa", "frequency") 
mk1 <- ggplot(df_kappa, aes(x=kappa, y=frequency)) + 
  geom_segment(aes(x=kappa, xend=kappa, y=0, yend=frequency/nout)) + 
  ylab("proportion")+ xlab(expression(kappa))
df_kappa <- data.frame(table(out_ppmx$kappangg[2,]))
colnames(df_kappa) <- c("kappa", "frequency")
mk2 <- ggplot(df_kappa, aes(x=kappa, y=frequency)) + 
  geom_segment(aes(x=kappa, xend=kappa, y=0, yend=frequency/nout)) + 
  ylab("proportion")+ xlab(expression(kappa))

ksp <- ggpubr::ggarrange(mk1, mk2, ms1, ms2, nrow=2, ncol = 2)


################################

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

#ccm <- matrix(0, nrow=10, ncol = 79)
ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:10){
  
  currfold <- (vectf[k]:(vectf[k+1]-1))
  X_train <- data.frame(X[-currfold,])
  Z_train <- data.frame(Z[-currfold,])
  Y_train <- data.frame(Y[-currfold,])
  
  X_test <- data.frame(X[currfold,])
  Z_test <- data.frame(Z[currfold,])
  Y_test <- data.frame(Y[currfold,])
  
  trtsgn_train <- trtsgn[-currfold]
  trtsgn_test <- trtsgn[currfold]
  
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)
  
  cls2 <- t(as.matrix(res0$label[[1]]))[,c(1:num_treat[1])]
  psm2 <- comp.psm(cls2)
  
  orp <- as.numeric(rownames(X[which(trtsgn == 1),]))
  id <- as.numeric(rownames(X_train[which(trtsgn_train == 1),]))
  
  #mc_vi2 <- minVI(psm2)
  #clu <- as.factor(mc_vi2$cl)
  co <- c()
  for(i in 1: num_treat[1]){
    co[i] <- which(orp == id[i])
  }
  #ccm[k,co] <- mc_vi2$cl
  ccm[co,co] <- ccm[co,co] + psm2
}

ccmf <- as.matrix(ccm/9)
#heatmap(ccmf)
colnames(ccmf) <- rownames(ccmf) <- orp
labels <- minVI(ccmf)$cl
labels <- minbinder(ccmf)$cl

reord <- c()
for(i in 1:max(labels)){
  reord <- c(reord, which(labels == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(data)
colnames(mC) <- c("Var1", "Var2", "Co-occurrence")
c1 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  theme_classic() +
  xlab("Patients") + ylab("Patients") #+ ggtitle("Heatmap of averaged co-occurence matrix")
c1
labels1 <- labels
#ccm <- matrix(0, nrow=10, ncol = 79)
ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:10){
  
  currfold <- (vectf[k]:(vectf[k+1]-1))
  X_train <- data.frame(X[-currfold,])
  Z_train <- data.frame(Z[-currfold,])
  Y_train <- data.frame(Y[-currfold,])
  
  X_test <- data.frame(X[currfold,])
  Z_test <- data.frame(Z[currfold,])
  Y_test <- data.frame(Y[currfold,])
  
  trtsgn_train <- trtsgn[-currfold]
  trtsgn_test <- trtsgn[currfold]
  
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)
  
  cls2 <- t(as.matrix(res0$label[[2]]))[,c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  
  orp <- as.numeric(rownames(X[which(trtsgn == 2),]))
  id <- as.numeric(rownames(X_train[which(trtsgn_train == 2),]))
  
  #mc_vi2 <- minVI(psm2)
  #clu <- as.factor(mc_vi2$cl)
  co <- c()
  for(i in 1: num_treat[2]){
    co[i] <- which(orp == id[i])
  }
  #ccm[k,co] <- mc_vi2$cl
  ccm[co,co] <- ccm[co,co] + psm2
}

ccmf <- as.matrix(ccm/9)
heatmap(ccmf)
colnames(ccmf) <- rownames(ccmf) <- orp
labels <- minVI(ccmf)$cl
labels <- minbinder(ccmf)$cl

reord <- c()
for(i in 1:max(labels)){
  reord <- c(reord, which(labels == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(data)
colnames(mC) <- c("Var1", "Var2", "Co-occurrence")
c2 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  theme_classic() +
  xlab("Patients") + ylab("Patients") #+ ggtitle("Heatmap of averaged co-occurence matrix")
c2

Y1 <- Y[which(trtsgn==1),]
rownames(Y1) <- rownames(matchRTComp)[which(trtsgn==1)]
apply(Y1[which(labels1 == 1),], 2, mean)
apply(Y1[which(labels1 == 2),], 2, mean)
apply(Y1[which(labels1 == 3),], 2, mean)

Y2 <- Y[which(trtsgn==2),]
rownames(Y2) <- rownames(matchRTComp)[which(trtsgn==2)]
apply(Y2[which(labels == 1),], 2, mean)
apply(Y2[which(labels == 2),], 2, mean)
apply(Y2[which(labels == 9),], 2, mean)
