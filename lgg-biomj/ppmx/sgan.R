rm(list=ls())
set.seed(121)
load(file = "output/lgg12aprs121.RData")
load("data/LGGdata.rda")
#name <- c("v01")

matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1

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

ccm <- matrix(0, nrow=10, ncol = 79)
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
#colnames(psm2) <- rownames(psm2) <- rownames(X_train[which(trtsgn_train == 2),])[reord2]
#heatmap(psm2)
mc_vi2 <- minVI(psm2)
orp <- as.numeric(rownames(X[which(trtsgn == 2),]))
id <- as.numeric(rownames(X_train[which(trtsgn_train == 2),]))
clu <- as.factor(mc_vi2$cl)
co <- c()
for(i in 1: num_treat[2]){
  co[i] <- which(orp == id[i])
}
ccm[k,co] <- mc_vi2$cl
}

av_psm <- comp.psm(ccm+1)
#minbinder.ext(av_psm)
avclu <- minVI(av_psm)$cl

X2 <- X[trtsgn==2,]
X2[avclu == 4,]
head(X2)


heatmap()
#reord2 <- c()
#for(i in 1:max(mc_vi2$cl)){
#  reord2 <- c(reord2, which(mc_vi2$cl == i))
#}
#tab <- table(mc_vi2$cl)
#print(tab[-1])
#rn <- rownames(X_train[which(trtsgn_train == 2),])[reord2]
#rn <- rn[-c(1:tab[1])]
#print(rn)



out_ppmx <- res0
# Co-occurence plot ----

data <- t(out_ppmx$label[[2]])
data <- data[,reord2]
coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
c2 <- ggplot(melt(coincidences), aes(Var1,Var2, fill=value/nout)) + geom_raster() +
  scale_fill_continuous(type = "viridis") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 2")

cp2 <- c2 + labs(fill = "Probability")
cp2

