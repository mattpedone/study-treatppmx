rm(list=ls())
set.seed(121)

library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

load("data/LGGdata.rda")
#name <- c("v01")

npc2 <- function(output, trtsgn, myoutot){
  K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- array(0, dim = c(3, 3, K))
  myctutSum <- NULL
  for (i in 1:K) {
    mycurdata <- output[, , i]
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
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
npat <- length(trtsgn)

K <- 2#repliche x convergenza
npat_pred <- 28

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

registerDoParallel(cores = K)#alloco solo core necessari

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

nout <- 2000

X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

myres0 <- foreach(k = 1:K) %dopar%
  {
    X_train <- data.frame(X[1:130,])
    Z_train <- data.frame(Z[1:130,])
    Y_train <- data.frame(Y[1:130,])
    
    X_test <- data.frame(X[131:158,])
    Z_test <- data.frame(Z[131:158,])
    Y_test <- data.frame(Y[131:158,])
    
    trtsgn_train <- trtsgn[1:130]
    trtsgn_test <- trtsgn[131:158]
    
    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_L0 <- diag(10, ncol(Y_train))
    modelpriors$hP0_nu0 <- ncol(Y_train) + 2; modelpriors$hP0_V0 <- diag(1, ncol(Y_train))
    
    #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 2000
    thinning <- 5
    
    nout <- (iterations-burnin)/thinning
    predAPT <- c()
    
    res0 <- tryCatch(expr = ppmxct(y = data.matrix(Y_train), X = data.frame(X_train), 
                                   Xpred = data.frame(X_test), Z = data.frame(Z_train), 
                                   Zpred = data.frame(Z_test), asstreat = trtsgn_train, #treatment,
                                   PPMx = 1, cohesion = 1, kappa = c(.01, 10, 5, 1), sigma = c(0.01, .5, 6),
                                   similarity = 2, consim = 2, similparam = vec_par, 
                                   calibration = 2, coardegree = 2, modelpriors, 
                                   update_hierarchy = T,
                                   hsp = T, iter = iterations, burn = burnin, thin = thinning, 
                                   mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, 
                                   nclu_init = 1), error = function(e){FALSE})
    return(res0)
  }

for(k in 1:K){
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[1:130]#simdata$trtsgn[[k]][1:124]
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
  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(myres[,, 1], myres[,, 2])
  
  nclust_all[k,] <- c(t(myclu))
  gof_all[k,] <- myfit
  
  #myprob <- simdata$prob[[k]]
}

myoutot <- as.numeric(matchRTComp[131:158,9])#simdata$yord[[k]][131:158,]
#NPC
PPMXCUT <- c()
for(k in 1:K){
  trtsgn_test <- trtsgn[131:158]#simdata$trtsgn[[k]][131:158]
  temp <- array(0, dim = c(28, 6, 1))
  temp[,,1] <- predAPT_all[, 4:9,k]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
}
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

PPMXRG <- c()
#ESM
#ho definito come respondent anche i partial responent
for(k in 1:K){
  mytab <- cbind(myass = predAPT_all[,3,k], rndass = trtsgn[131:158], resp = as.numeric(myoutot>=2))
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
  if(length(row.names(table1)) == 2){
    crt1 <- table1[2,1]/sum(table1[,1])
  }
  
  if(length(table2) == 4){
    crt2 <- table2[2,2]/sum(table2[,2])
  }
  if(length(table2) < 4){
    crt2 <- as.numeric(row.names(table2))
  }
  if(length(row.names(table2)) == 2){
    crt2 <- table2[2,1]/sum(table2[,1])
  }
  

  ### summary meaures
  PPMXRG[k] <- crt1*p1 + crt2*p2 - sum(as.numeric(myoutot>=2))/npat_pred
}

ESM <- c(round(mean(PPMXRG), 4), round(sd(PPMXRG), 4))

#FIT
WAIC <- c(mean(gof_all[,1]), sd(gof_all[,1]))
lpml <- c(mean(gof_all[,2]), sd(gof_all[,2]))

#results
resPPMX <- rbind(NPC, ESM, WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")
resPPMX

cluPPMX <- nclust_all[,-c(3,4)]
clu <- apply(cluPPMX, 2, mean)
clu <- rbind(clu, apply(cluPPMX, 2, sd))
colnames(clu) <- c("avg # trt 1", "avg # trt 2", "VI trt 1", "VI trt 2")
#clu

NPC; ESM; PPMXCUT[1]; PPMXRG[1]
#save(resPPMX, file=paste0(mypath, sc, "res.RData"))
#save(cluPPMX, file=paste0(mypath, sc, "clu.RData"))
#save(PPMXCT, file=paste0(mypath, sc, "mot.RData"))
#save(PPMXpp, file=paste0(mypath, sc, "mtug.RData"))
#save(PPMXCUT, file=paste0(mypath, sc, "npc.RData"))
#}

out_ppmx <- res0

mc <- apply(out_ppmx$nclu, 1, mean)
trt <- trtsgn[1:130]#trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
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

data <- t(out_ppmx$label[[1]])
data <- data[,reord]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(coincidences)
c1 <- ggplot(mC, aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  xlab("patients") + ylab("patients") + ggtitle("Treatment 1")

cp1 <- c1 + labs(fill = "Correlation")

data <- t(out_ppmx$label[[2]])
data <- data[,reord2]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
c2 <- ggplot(melt(coincidences), aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  xlab("patients") + ylab("patients") + ggtitle("Treatment 2")

cp2 <- c2 + labs(fill = "Correlation")

#corrp_vi <- ggpubr::ggarrange(cp1, cp2, nrow=1, ncol = 2, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
#corrp_vi
#ggsave(corrp, device = "pdf", path = "figs", filename = "corr_plot.pdf")

data <- t(out_ppmx$label[[1]])
data <- data[,reordb]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(coincidences)
c1 <- ggplot(mC, aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  xlab("patients") + ylab("patients") + ggtitle("Treatment 1")

cp1b <- c1 + labs(fill = "Correlation")

data <- t(out_ppmx$label[[2]])
data <- data[,reord2b]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
c2 <- ggplot(melt(coincidences), aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  xlab("patients") + ylab("patients") + ggtitle("Treatment 2")

cp2b <- c2 + labs(fill = "Correlation")

#corrp_b <- ggpubr::ggarrange(cp1, cp2, nrow=1, ncol = 2, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
corr <- ggpubr::ggarrange(cp1, cp2, cp1b, cp2b, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
corr

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

# Posterior frequency for (\kappa, \sigma) ----

par(mfrow=c(2,1))
barplot(table(out_ppmx$sigmangg[1,]), breaks = 10, xlab = expression(sigma),
        main = expression(paste(sigma, " Marginal - Treatment 1")))
barplot(table(out_ppmx$sigmangg[2,]), breaks = 10, xlab = expression(sigma),
        main = expression(paste(sigma, " Marginal - Treatment 2")))

#dev.print(pdf, "figs/marg_sigma.pdf") 

#plot(out_ppmx$sigmangg[1,], type ="l")
#plot(out_ppmx$sigmangg[2,], type ="l")
#table(out_ppmx$sigmangg[1,])
#table(out_ppmx$sigmangg[2,])

barplot(table(out_ppmx$kappangg[1,]), breaks = 10, xlab = expression(kappa),
        main = expression(paste(kappa, " Marginal - Treatment 1")))
barplot(table(out_ppmx$kappangg[2,]), breaks = 10, xlab = expression(kappa),
        main = expression(paste(kappa, " Marginal - Treatment 2")))

apply(out_ppmx$beta, c(1, 2), mean)

myplot <- function(df){
  colnames(df) <- c("beta")
  df <- cbind(Index = as.numeric(row.names(df)), df)
  df <- reshape2::melt(df, id.vars="Index")
  colnames(df) <- c("Index", "beta", "value")
  tp <- ggplot2::ggplot(df, aes(x = Index, y = value)) + 
    geom_line() + theme_classic() + 
    xlab("Iterations") + ylab("beta") 
  
  den <- ggplot(df, aes(x=value)) + geom_density() +  
    geom_vline(aes(xintercept = quantile(value, probs = .025)), color="blue", 
               linetype="dashed", size=.25) + 
    geom_vline(aes(xintercept = quantile(value, probs = .975)), color="blue", 
               linetype="dashed", size=.25) + 
    geom_vline(aes(xintercept = 0), color="red", linetype="dashed", size=.25) + 
    xlab("beta")
  
  p <- ggpubr::ggarrange(tp, den, nrow = 1, ncol = 2)
  return (p)
}

b11 <- myplot(data.frame(out_ppmx$beta[1,1,]))
b12 <- myplot(data.frame(out_ppmx$beta[1,2,]))
b13 <- myplot(data.frame(out_ppmx$beta[1,3,]))
b21 <- myplot(data.frame(out_ppmx$beta[2,1,]))
b22 <- myplot(data.frame(out_ppmx$beta[2,2,]))
b23 <- myplot(data.frame(out_ppmx$beta[2,3,]))

progn <- ggpubr::ggarrange(b11, b12, b13, b21, b22, b23, nrow=2, ncol = 3, 
                           common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
progn

#posterior predictive probabilities ----
A0 <- apply(out_ppmx$ypred, c(1,2,3), mean, na.rm=TRUE);#A0
c(as.numeric(c(A0[,,1]%*%wk)<c(A0[,,2]%*%wk))+1)

#A0 <- c(apply(out_ppmx$pipred, c(1,2,3), median, na.rm=TRUE))#, mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)

#posterior distribution of predictive utility
ns <- dim(out_ppmx$pipred)[4]
npp <- dim(out_ppmx$pipred)[1]
for(pt in 1:npp){
  dpu <- matrix(0, ns, 2)
  for(i in 1:ns){
    dpu[i,] <- apply(out_ppmx$pipred[pt,,,i]*wk, 2, sum)
  }
  
  df <- data.frame(dpu)
  colnames(df) <- c("treatment 1", "treatment 2")
  df <- cbind(Index = as.numeric(row.names(df)), df)
  df <- reshape2::melt(df, id.vars="Index")
  colnames(df) <- c("Index", "Treatment", "value")
  pt1 <- ggplot2::ggplot(df, aes(x = value, col = Treatment)) + 
    geom_histogram(alpha = 0.5, position = "identity") + theme_classic() + 
    xlab("Predicted Utility") + ylab("") + ggtitle(paste0("Patient", pt))
  assign(paste("pt", pt, sep=""),pt1)
}

ut <- ggpubr::ggarrange(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11,
                        pt12, pt13, pt14, pt15, pt18, nrow=4, ncol = 4, 
                        common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
ut
