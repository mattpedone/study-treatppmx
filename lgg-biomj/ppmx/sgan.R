rm(list=ls())
library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(reshape2)
library(ggridges)
library(viridis)
library(tidyverse)

library(plotly)
library(dplyr)
set.seed(121)
load(file = "output/lgg12aprs121.RData")
load("data/LGGdata.rda")
#name <- c("v01")

matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
mycn <- c(colnames(matchRTComp)[1:10], "ACVRL1", "Src", "HSP70", 
          "HER2", "PAI-1", "Smad1", "FOXO3a", 
          "PRAS40", "Cyclin", "SF2", "Bad" , 
          "Lck", "Caspase-7_cleavedD198", "Paxillin", "MYH11", 
          "14-3-3_epsilon", "Akt", "Caveolin-1", "Rab25", 
          "YAP", "RBM15", "Claudin-7", "ER-alpha", 
          "C-Raf", "CD31", "Ku80", "Bcl-2", "GSK3-alpha-beta")
colnames(matchRTComp) <- mycn

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
g5 <- orp[which(labels==1)]

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
#ggsave(c1, device = "pdf", path = "figs", filename = "avg_coocc.pdf")

#PREDITTIVE
pred_cov <- matchRTComp[, c(16:38)]
pred_cov_g1 <- pred_cov[as.character(g1),]
pred_cov_g2 <- pred_cov[as.character(g2),]
pred_cov_g3 <- pred_cov[as.character(g3),]
pred_cov_g4 <- pred_cov[as.character(g4),]
pred_cov_g5 <- pred_cov[as.character(g5),]

#par(mfrow = c(2,2))
#
#for(ind in 16:38){
#plot(density(matchRTComp[matchRTComp$newTRT==1, ind]), 
#     xlab = colnames(matchRTComp)[ind],
#     main = " ")#paste0("Empirical dens ", colnames(matchRTComp)[ind]))
#rug(jitter(pred_cov_g1[,ind-15]),col="blue",lwd=2)
#rug(jitter(pred_cov_g2[,ind-15]),col="red",lwd=2)
#rug(jitter(pred_cov_g3[,ind-15]),col="green",lwd=2)
#rug(jitter(pred_cov_g4[,ind-15]),col="magenta",lwd=2)
##rug(jitter(pred_cov_g5[,ind-15]),col="orange",lwd=2)
#}
#
#cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#fino a 38
df_pred1 <- data.frame(matchRTComp[matchRTComp$newTRT==1, c(16:24)])
df_pred1$Group <- as.character(labels)
preddens <- df_pred1 %>%
  pivot_longer(!Group, names_to = "variable", values_to = "Protein Expression") %>%
  ggplot(aes(x = `Protein Expression`)) + geom_density() + 
  #geom_point(aes(x =GSK3.alpha.beta, y= 0), col = labels, shape=labels, size=2.5) 
  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc"))+
  facet_wrap(~variable, scales = "free", ncol = 3)
  #ggtitle("Predictive covariates densities") #+ 
ggsave(preddens, device = "pdf", path = "figs", filename = "pred_dens1.pdf")

df_pred2 <- data.frame(matchRTComp[matchRTComp$newTRT==1, c(25:33)])
df_pred2$Group <- as.character(labels)
preddens <- df_pred2 %>%
  pivot_longer(!Group, names_to = "variable", values_to = "Protein Expression") %>%
  ggplot(aes(x = `Protein Expression`)) + geom_density() + 
  #geom_point(aes(x =GSK3.alpha.beta, y= 0), col = labels, shape=labels, size=2.5) 
  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc"))+
  facet_wrap(~variable, scales = "free", ncol = 3)
ggsave(preddens, device = "pdf", path = "figs", filename = "pred_dens2.pdf")

df_pred2 <- data.frame(matchRTComp[matchRTComp$newTRT==1, c(34:38)])
df_pred2$Group <- as.character(labels)
preddens <- df_pred2 %>%
  pivot_longer(!Group, names_to = "variable", values_to = "Protein Expression") %>%
  ggplot(aes(x = `Protein Expression`)) + geom_density() + 
  #geom_point(aes(x =GSK3.alpha.beta, y= 0), col = labels, shape=labels, size=2.5) 
  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc"))+
  facet_wrap(~variable, scales = "free", ncol = 3)
ggsave(preddens, device = "pdf", path = "figs", filename = "pred_dens3.pdf")
  
df_poster <- data.frame(matchRTComp[matchRTComp$newTRT==1, c(22, 31, 36)])
df_poster$Group <- as.character(labels)
preddens <- df_poster %>%
  pivot_longer(!Group, names_to = "variable", values_to = "Protein Expression") %>%
  ggplot(aes(x = `Protein Expression`)) + geom_density() + 
  #geom_point(aes(x =GSK3.alpha.beta, y= 0), col = labels, shape=labels, size=2.5) 
  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc"))+
  facet_wrap(~variable, scales = "free", ncol = 1) + 
  theme_minimal() 

poster_grid <- cowplot::plot_grid(c1, preddens, align = "h", ncol = 2, rel_widths = c(5/8, 3/8))
ggsave(poster_grid, device = "pdf", path = "../../conferences/isba22/", filename = "pred_poster.pdf")
#rg5 <- as.character(myc[as.character(g5),2])
#lab_g5 <- as.character(lab[as.character(g5),])

#PROGNOSTICHE
prog_cov <- matchRTComp[, c(11, 13)]
prog_cov_g1 <- prog_cov[as.character(g1),]
prog_cov_g2 <- prog_cov[as.character(g2),]
prog_cov_g3 <- prog_cov[as.character(g3),]
prog_cov_g4 <- prog_cov[as.character(g4),]
prog_cov_g5 <- prog_cov[as.character(g5),]

#par(mfrow = c(1, 2))
#plot(density(matchRTComp[matchRTComp$newTRT==1, c(11)]), 
#     xlab = colnames(matchRTComp)[11],
#     main = paste0("Empirical dens ", colnames(matchRTComp)[11]))
#rug(jitter(prog_cov_g1$`ACVRL1-R-C`),col="blue",lwd=2)
#rug(jitter(prog_cov_g2$`ACVRL1-R-C`),col="red",lwd=2)
#rug(jitter(prog_cov_g3$`ACVRL1-R-C`),col="green",lwd=2)
#rug(jitter(prog_cov_g4$`ACVRL1-R-C`),col="magenta",lwd=2)
##rug(jitter(prog_cov_g5$`ACVRL1-R-C`),col="orange",lwd=2)

#plot(density(matchRTComp[matchRTComp$newTRT==1, c(13)]), 
#     xlab = colnames(matchRTComp)[13],
#     main = paste0("Empirical dens ", colnames(matchRTComp)[13]))
#rug(jitter(prog_cov_g1$`HSP70-R-C`),col="blue",lwd=2)
#rug(jitter(prog_cov_g2$`HSP70-R-C`),col="red",lwd=2)
#rug(jitter(prog_cov_g3$`HSP70-R-C`),col="green",lwd=2)
#rug(jitter(prog_cov_g4$`HSP70-R-C`),col="magenta",lwd=2)
##rug(jitter(prog_cov_g5$`HSP70-R-C`),col="orange",lwd=2)

df_prog <- data.frame(matchRTComp[matchRTComp$newTRT==1, c(11, 13)])
df_prog$Group <- as.character(labels)
progdens <- df_prog %>%
  pivot_longer(!Group, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) + geom_density() + 
  #geom_point(aes(x =GSK3.alpha.beta, y= 0), col = labels, shape=labels, size=2.5) 
  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc"))+
  facet_wrap(~variable, scales = "free", ncol = 2) + 
  ggtitle("Prognostic covariates densities") #+ 
ggsave(progdens, device = "pdf", path = "figs", filename = "prog_dens.pdf")

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
prob_5 <- exp(as.matrix(prog_cov_g5)%*%beta)

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
lab_g5 <- as.character(lab[as.character(g5),])

df1 <- data.frame(prob_p = round(prob_1/rowSums(prob_1), 4), 
                  prob_pp = round(predprob[as.character(g1),], 4), response = lab_g1)
df2 <- data.frame(prob_p = round(prob_2/rowSums(prob_2), 4), 
                  prob_pp = round(predprob[as.character(g2),], 4), response = lab_g2)
df3 <- data.frame(prob_p = round(prob_3/rowSums(prob_3), 4), 
                  prob_pp = round(predprob[as.character(g3),], 4), response = lab_g3)
df4 <- data.frame(prob_p = round(prob_4/rowSums(prob_4), 4), 
                  prob_pp = round(predprob[as.character(g4),], 4), response = lab_g4)
df5 <- data.frame(prob_p = round(prob_5/rowSums(prob_5), 4), 
                  prob_pp = round(predprob[as.character(g5),], 4), response = lab_g5)

#gruppo 1 
df1 <- cbind(df1, var_cc=c(0.3990, -0.1455, -0.4622))
#gruppo 2 
df2 <- cbind(df2, var_cc=c(df2[,5]-df2[,2])/df2[,2])
#gruppo 3 
df3 <- cbind(df3, var_cc=c(0.3278, -0.0694, -0.0691, 0.0031, 0.3513, -0.0255))
#gruppo 4
df4 <- cbind(df4, var_cc=c(0.1553, 0.2041, 0.1776, 0.2891, -0.3054, 0.4334, 
                           0.0381, 0.1729, 0.0139, 0.2861, 0.1735, -0.1075, 
                           -0.0265))

#SIMILARITY MATRIX - DENSITIES
#sim_mat <- matrix(NA, 79, 79)
##pred_cov2 <- scale(pred_cov[which(trtsgn == 2),])
#pred_cov2 <- pred_cov[which(trtsgn == 2),]
#for(i in 1:78){
#  for(j in (i+1):79){
#    sim <- 0
#    for(p in 1:23){
#      somma <- pred_cov2[i, p] + pred_cov2[j, p]
#      a <- treatppmx::gsimconNNIG(m0 = 0, k0 = 1, nu0 = 2.0, s20 = 1.0, 
#                                  sumx = somma, sumx2 = somma*somma, 
#                                  n = 2, DD = 1, logout = 0)
#      sim <- sim + a
#    }
#    sim_mat[i, j] <- sim
#  }
#}
#
#colnames(sim_mat) <- rownames(sim_mat)<- as.character(orp)
#mat <- melt(sim_mat)
#mat <- subset(mat, value != 0.0)
#con <- rownames_to_column(data.frame(Group = t(t(labels))))
#con$Group <- as.factor(con$Group)
#mat$Group <- NA
#for(i in 1:dim(mat)[1]) {
#  val <- con$Group[which(as.numeric(con$rowname) == mat$Var1[i])]
#  mat$Group[i] <- val
#}
#
#df <- data.frame(mat)
#df$Group <- as.factor(df$Group)
#colnames(df)[3] <- c("Similarity")
#dens <- df %>%
#  ggplot(aes(x = Similarity)) + geom_density() + theme(legend.position="none") + 
#  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc")) 
#simdens <- dens
#
#df <- data.frame(mat)
#df$Group <- as.factor(df$Group)
#colnames(df)[3] <- c("Similarity")
#dens <- df %>%
#  ggplot(aes(x = Similarity)) + stat_ecdf() + ylab("cumulative probability") + #theme(legend.position="none") +
#  geom_rug(mapping = aes(color = Group), sides = "b", size = 1, length = unit(0.05, "npc")) 
#simecdf <- dens
#
#sim_grid <- cowplot::plot_grid(simdens, simecdf, align = "h", ncol = 2, 
#                               labels = c("A", "B"), rel_widths = c(15/32, 17/32))
#sim_grid

#SIMILARITY MATRIX
#pred_cov2 <- scale(pred_cov[which(trtsgn == 2),])
pred_cov2 <- pred_cov[which(trtsgn == 2),]

sim_mat <- matrix(NA, dim(pred_cov2)[1], dim(pred_cov2)[1])

for(i in 1:dim(pred_cov2)[1]){
  for(j in 1:dim(pred_cov2)[1]){
    sim <- 1
    for(p in 1:dim(pred_cov2)[2]){
      somma <- pred_cov2[i, p] + pred_cov2[j, p]
      sim <- sim + treatppmx::gsimconNNIG(m0 = 0, k0 = 1, nu0 = 2.0, s20 = 1.0, 
                                          sumx = somma, sumx2 = somma*somma, 
                                          n = 2, DD = 2, logout = 1)
      sim <- sim+(1/sqrt(dim(pred_cov2)[2]))
    }
    sim_mat[i, j] <- sim
  }
}
reord <- c()
for(i in 1:max(labels)){
  reord <- c(reord, which(labels == i))
}

# Co-occurence plot ----
data1 <- sim_mat[,reord]
sim_mat <- data1[reord,]
#colnames(data) <- rownames(data) <- NULL
#sim_mat[lower.tri(sim_mat)] <- NA
#diag(sim_mat) <- NA
mc2 <- melt(sim_mat)
colnames(mc2) <- c("Var1", "Var2", "Similarity")
simp <- ggplot(mc2, aes(Var1, Var2, fill=`Similarity`)) + geom_tile() + #color = "white"
  #filter(row_number() >= which(Var1 == Var2)) +
  #geom_tile(aes(fill = value), color='white') +
  scale_fill_continuous(type = "viridis") + 
  xlab("Patients") + ylab("Patients") + 
  theme_classic() +
  theme(panel.grid = element_blank()) + #+ ggtitle("Heatmap of averaged co-occurence matrix")
  geom_rect(mapping = aes(xmin = 0.5, xmax = 53.5, ymin = 0.5, ymax = 53.5),
          fill = NA, col = "black") +
  geom_rect(mapping = aes(xmin = 53.5, xmax = 66.5, ymin = 53.5, ymax = 66.5),
            fill = NA, col = "black") +
  geom_rect(mapping = aes(xmin = 66.5, xmax = 69.5, ymin = 66.5, ymax = 69.5),
          fill = NA, col = "black") + 
  geom_rect(mapping = aes(xmin = 69.5, xmax = 73.5, ymin = 69.5, ymax = 73.5),
          fill = NA, col = "black") + 
  geom_rect(mapping = aes(xmin = 73.5, xmax = 79.5, ymin = 73.5, ymax = 79.5),
            fill = NA, col = "black")

ggsave(simp, device = "pdf", path = "figs", filename = "sim_plot.pdf")

X[which(labels==3),]


####################################
labels <- minVI(ccmf)$cl

g1 <- orp[which(labels==3)]
g2 <- orp[which(labels==4)]
g3 <- orp[which(labels==5)]
g4 <- orp[which(labels==2)]
g0 <- orp[which(labels==1)]

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
prob_5 <- exp(as.matrix(prog_cov_g5)%*%beta)

predproball <- c()
for(k in 1:10){
  res0 <- myres0[[k]]
  restmp <- apply(res0$ypred, c(1, 2, 3), mean)
  trttmp <- as.numeric(restmp[,,1]%*%wk<restmp[,,2]%*%wk+1)+1
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
