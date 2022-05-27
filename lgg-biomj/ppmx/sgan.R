rm(list=ls())
library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)
set.seed(121)
load(file = "output/lgg12aprs121.RData")
load("data/LGGdata.rda")
#name <- c("v01")

matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
wk <- c(0, 40, 100)
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)
#load("/home/matt/Dropbox/PHD/study-treatppmx/output/lgg_analysis_24mar.RData")
nout <- 1600
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

npat <- length(trtsgn)

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
colnames(ccmf) <- rownames(ccmf) <- orp
#hmccm <- heatmap(ccmf)
#
#hmccmc <- ccmf
#hmccmc1 <- hmccmc[,hmccm$rowInd]
#hmccmc2 <- hmccmc1[hmccm$colInd,]
#
#heatmap(hmccmc2, Colv = NA, labRow = NA, main = "Heatmap of averaged co-occurence matrix")
#
#hc <- hclust(dist(ccmf), method = "ave")
##cutree(hc, h=1.5)
#plot(hc, xlab = " ", main = "Heatmap dendogram")
#abline(h=1.5, col = "red")
#
#g1 <- orp[which(cutree(hc, h=1.5)==2)]
#g2 <- orp[which(cutree(hc, h=1.5)==3)]
#g3 <- orp[which(cutree(hc, h=1.5)==4)]
#g4 <- orp[which(cutree(hc, h=1.5)==5)]

labels <- minVI(ccmf)$cl

g1 <- orp[which(labels==3)]
g2 <- orp[which(labels==4)]
g3 <- orp[which(labels==5)]
g4 <- orp[which(labels==2)]

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
c1 <- ggplot(mC, aes(Var1,Var2, fill=value)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  xlab("Patients") + ylab("Patients") + ggtitle("Heatmap of averaged co-occurence matrix")
c1

#PREDITTIVE
pred_cov <- matchRTComp[, c(16:38)]
pred_cov_g1 <- pred_cov[as.character(g1),]
pred_cov_g2 <- pred_cov[as.character(g2),]
pred_cov_g3 <- pred_cov[as.character(g3),]
pred_cov_g4 <- pred_cov[as.character(g4),]
#pred_cov_g5 <- pred_cov[as.character(g5),]

par(mfrow = c(3,4))

for(ind in 16:38){
plot(density(matchRTComp[matchRTComp$newTRT==1, ind]), 
     xlab = colnames(matchRTComp)[ind],
     main = " ")#paste0("Empirical dens ", colnames(matchRTComp)[ind]))
rug(jitter(pred_cov_g1[,ind-15]),col="blue",lwd=2)
rug(jitter(pred_cov_g2[,ind-15]),col="red",lwd=2)
rug(jitter(pred_cov_g3[,ind-15]),col="green",lwd=2)
rug(jitter(pred_cov_g4[,ind-15]),col="magenta",lwd=2)
#rug(jitter(pred_cov_g5[,ind-15]),col="orange",lwd=2)
}



#rg5 <- as.character(myc[as.character(g5),2])
#lab_g5 <- as.character(lab[as.character(g5),])

#PROGNOSTICHE
prog_cov <- matchRTComp[, c(11, 13)]
prog_cov_g1 <- prog_cov[as.character(g1),]
prog_cov_g2 <- prog_cov[as.character(g2),]
prog_cov_g3 <- prog_cov[as.character(g3),]
prog_cov_g4 <- prog_cov[as.character(g4),]
#prog_cov_g5 <- prog_cov[as.character(g5),]

par(mfrow = c(1, 2))
plot(density(matchRTComp[matchRTComp$newTRT==1, c(11)]), 
     xlab = colnames(matchRTComp)[11],
     main = paste0("Empirical dens ", colnames(matchRTComp)[11]))
rug(jitter(prog_cov_g1$`ACVRL1-R-C`),col="blue",lwd=2)
rug(jitter(prog_cov_g2$`ACVRL1-R-C`),col="red",lwd=2)
rug(jitter(prog_cov_g3$`ACVRL1-R-C`),col="green",lwd=2)
rug(jitter(prog_cov_g4$`ACVRL1-R-C`),col="magenta",lwd=2)
#rug(jitter(prog_cov_g5$`ACVRL1-R-C`),col="orange",lwd=2)

plot(density(matchRTComp[matchRTComp$newTRT==1, c(13)]), 
     xlab = colnames(matchRTComp)[13],
     main = paste0("Empirical dens ", colnames(matchRTComp)[13]))
rug(jitter(prog_cov_g1$`HSP70-R-C`),col="blue",lwd=2)
rug(jitter(prog_cov_g2$`HSP70-R-C`),col="red",lwd=2)
rug(jitter(prog_cov_g3$`HSP70-R-C`),col="green",lwd=2)
rug(jitter(prog_cov_g4$`HSP70-R-C`),col="magenta",lwd=2)
#rug(jitter(prog_cov_g5$`HSP70-R-C`),col="orange",lwd=2)

#baseline probabilities & predicted probabilities
beta <- array(0, dim=c(2, 3, 10))
for(k in 1:10){
  res0 <- myres0[[k]]
  beta[,, k] <- apply(res0$beta, c(1, 2), mean)
}

beta <- apply(beta, c(1, 2), mean)

prob_1 <- exp(as.matrix(prog_cov_g1)%*%beta)
prob_2 <- exp(as.matrix(prog_cov_g2)%*%beta)
prob_3 <- exp(as.matrix(prog_cov_g3)%*%beta)
prob_4 <- exp(as.matrix(prog_cov_g4)%*%beta)
#prob_5 <- exp(as.matrix(prog_cov_g5)%*%beta)

predproball <- c()
for(k in 1:10){
  res0 <- myres0[[k]]
  restmp <- apply(res0$ypred, c(1, 2, 3), mean)
  trttmp <- as.numeric(restmp[,,1]%*%wk<restmp[,,2]%*%wk)+1
  predprob <- matrix(0, nrow(restmp), 3)
  for(i in 1:nrow(restmp)){
    predprob[i,] <- restmp[i, , trttmp[i]]
  }
  predproball <- rbind(predproball, predprob)
}

predprob <- predproball[which(trtsgn == 2),]
rownames(predprob) <- as.character(orp)

#relabel response
myc <- matchRTComp[, 2:10]
lab <- matchRTComp[, 3]
lab <- data.frame(dplyr::recode(lab, "Stable Disease" = "PS", "Partial Remission/Response" = "PS", 
                                "Complete Remission/Response" = "CR", "Progressive Disease" = "PD"))
rownames(lab) <- rownames(myc)
lab_g1 <- as.character(lab[as.character(g1),])
lab_g2 <- as.character(lab[as.character(g2),])
lab_g3 <- as.character(lab[as.character(g3),])
lab_g4 <- as.character(lab[as.character(g4),])

df1 <- data.frame(prob_p = round(prob_1/rowSums(prob_1), 4), 
                  prob_pp = round(predprob[as.character(g1),], 4), response = lab_g1)
df2 <- data.frame(prob_p = round(prob_2/rowSums(prob_2), 4), 
                  prob_pp = round(predprob[as.character(g2),], 4), response = lab_g2)
df3 <- data.frame(prob_p = round(prob_3/rowSums(prob_3), 4), 
                  prob_pp = round(predprob[as.character(g3),], 4), response = lab_g3)
df4 <- data.frame(prob_p = round(prob_4/rowSums(prob_4), 4), 
                  prob_pp = round(predprob[as.character(g4),], 4), response = lab_g4)
#gruppo 1 
df1
#gruppo 2 
df2 
#gruppo 3 
df3
#gruppo 4
df4

#SIMILARITY MATRIX
sim_mat <- matrix(0, 79, 79)
pred_cov2 <- pred_cov[which(trtsgn == 2),]
for(i in 1:78){
  for(j in (i+1):79){
    sim <- 0
    for(p in 1:23){
      somma <- pred_cov2[i, p] + pred_cov2[j, p]
      sim <- sim + treatppmx::gsimconNNIG(m0 = 0, k0 = 1, nu0 = 2.0, s20 = 1.0, 
                                          sumx = somma, sumx2 = somma*somma, 
                                          n = 2, DD = 1, logout = 0)
    }
    sim_mat[i, j] <- sim
  }
}
colnames(sim_mat) <- rownames(sim_mat)<- as.character(orp)

d1 <- c(sim_mat[c(62, 7), 72], sim_mat[c(7), 62])
d2 <- c(sim_mat[c(10, 44, 73), 79], sim_mat[c(10, 44), 73], sim_mat[c(10), 44])
d3 <- c(sim_mat[c(20, 24, 60, 61, 77), 78], sim_mat[c(20, 24, 60, 61), 77], 
        sim_mat[c(20, 24, 60), 61], sim_mat[c(20, 24), 60], sim_mat[20, 24])
d4 <- c(sim_mat[c(5, 14, 16, 17, 38, 40, 41, 46, 55, 63, 67), 74], 
        sim_mat[c(5, 14, 16, 17, 38, 40, 41, 46, 55, 63), 67], 
        sim_mat[c(5, 14, 16, 17, 38, 40, 41, 46, 55), 63], 
        sim_mat[c(5, 14, 16, 17, 38, 40, 41, 46), 55], 
        sim_mat[c(5, 14, 16, 17, 38, 40, 41), 46], 
        sim_mat[c(5, 14, 16, 17, 38, 40), 41], 
        sim_mat[c(5, 14, 16, 17, 38), 40], sim_mat[c(5, 14, 16, 17), 38],  
        sim_mat[c(5, 14, 16), 17], sim_mat[c(5, 14), 16], sim_mat[5, 14])
plot(density(sim_mat[upper.tri(sim_mat, diag=FALSE)]), main = "emp dist")
#plot(ecdf(sim_mat[upper.tri(sim_mat, diag=FALSE)]), main = "ecdf similarity distances")
rug(jitter(d4),col="magenta",lwd=2)
rug(jitter(d1),col="blue",lwd=2)
rug(jitter(d2),col="red",lwd=2)
rug(jitter(d3),col="green",lwd=2)

plot(ecdf(sim_mat[upper.tri(sim_mat, diag=FALSE)]), main = "ecdf similarity distances")
rug(jitter(d4),col="magenta",lwd=2)
rug(jitter(d1),col="blue",lwd=2)
rug(jitter(d2),col="red",lwd=2)
rug(jitter(d3),col="green",lwd=2)

