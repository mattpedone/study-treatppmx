library("Rcpp")
Rcpp::sourceCpp('calcola_dist.cpp')

set.seed(123)
X <- matrix(rnorm(40),ncol=4)
X[,3] <- rbinom(5,p=0.4,1)
X[,4] <- rbinom(5,p=0.4,1)

head(X)

dist_r(X,2,2)
## [1] 42

centroide(X,2,2)




alldata <- read.table("~/Vecchio_Desktop/ESEMPIO_AVIS_DATIPPMX/build/dtasim.txt")
head(alldata)

Indici <- read.table("~/Vecchio_Desktop/ESEMPIO_AVIS_DATIPPMX/build/indici.txt")

Indici <- as.matrix(Indici)
out <- vector(length=100)
for(idx in 1:100){
data<- alldata[Indici[idx,],]
n = dim(data)[1]
p = dim(data)[2]
X = as.matrix(data[, 2:p])
out[idx] <- calcola_D_norm(X,1,2);
}

hist(out)

XX <- as.matrix(alldata[,c(2,3,4)])

distanze=dist_r(XX,1,2)
#plot(distanze$dist)
calcola_D(XX,1,2);
calcola_D_norm(XX,1,2);
calcola_D(XX[-c(1,2),],1,2);



Nj <- 2:200


Dmeannj <- matrix(nrow=200,ncol=3) 
Dmeannj[1] <- 0

G <- 5000
for(j in Nj){
Dmean <- vector(length=G)
for(g in 1:G){
idx <- sample(1:1000,size=j,replace=F)
Dmean[g] <- calcola_D_norm(XX[idx,],1,2);
}

Dmeannj[j,]= quantile(Dmean,prob=c(0.025,0.5,0.975))

}



apply(Dmeannj,2,diff)

#####
Dmeannjold <- matrix(nrow=200,ncol=3)
Dmeannj[1] <- 0


G <- 5000
for(j in Nj){
Dmean <- vector(length=G)
for(g in 1:G){
idx <- sample(1:1000,size=j,replace=F)
Dmean[g] <- calcola_D(XX[idx,],1,2);
}

Dmeannjold[j,]=quantile(Dmean,prob=c(0.025,0.5,0.975))
}

#x11()

par(mfrow=c(1,2))
matplot(1:200,100*Dmeannj,type="l")

matplot(1:200,100*Dmeannjold,type="l")



 
 x11()
 
 par(mfrow=c(1,2))
 matplot(100*apply(Dmeannj,2,diff),type="l")
 abline(h=0)
 matplot(0.037*apply(Dmeannjold,2,diff),type="l")
 abline(h=0)
# 
# hist(Dmean)
# mean(Dmean)
# var(Dmean)
# min(Dmean)
# max(Dmean)
# 
# 
# plot(Dmeannj/mean(diff(Dmeannj)))
# 
# Xu <- unique(XX)
# 
# Nu <- dim(Xu)[1]
# nj <- rep(0,Nu)
# JJ <- list()
# for(j in 1:Nu){
# 	JJ[[j]]=1
# 	for(i in 1:1000){
# 		nj[j] =nj[j]+all(Xu[j,]==XX[i,])
# 		if(all(Xu[j,]==XX[i,]))
# 		{JJ[[j]]=c(JJ[[j]],i)}
# 	}
# 	JJ[[j]]=JJ[[j]][-1]
# }
# 
# 
# JJ
# 
# out <- vector(length=G)
# for(g in 1:G){
# iii <- sample(1:1000,1)
# out[g] <- calcola_D(rbind(XX[JJ[[1]],],XX[iii,]),1,2)
# }
# 
# 
# hist(out)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# Nj <- 65
# idx <- sample(1:1000,size=Nj,replace=F)
# 
# Xe <- XX[idx,]
# d <- calcola_D(Xe,1,2);
# d
# 
# out <- vector(length=G)
# for(g in 1:G){
# iii <- sample(1:1000,1)
# out[g] <- calcola_D(rbind(Xe,XX[iii,]),1,2)
# }
# 
# 
# hist(out)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# alldataj <- read.table("~/Desktop/ESEMPIO_AVIS_DATIPPMX/build/dtasimj.txt")
# 
# XXj <- as.matrix(alldataj[,c(2,3,4)])
# head(XXj)
# distanzej=dist_r(XXj,1,2)
# 
# length(unique(distanzej$dist))
# 
# plot(distanzej$dist-distanze$dist)
# 
# 
