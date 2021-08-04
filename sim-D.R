library("Rcpp")
Rcpp::sourceCpp("src/calcola_dist.cpp")
load("/home/matt/Dropbox/PHD/study-treatppmx/data/SimuOutsce2.rda")

XX <- t(mydata)[,c(1:10)]
n <- nrow(XX)
P <- ncol(XX)

Nj <- 10

Dmeannj <- matrix(nrow = Nj, ncol = 3) 
Dmeannj[1,] <- 0
G <- 5000
#G <- 50
for(j in 2:Nj){
  Dmean <- vector(length=G)
  for(g in 1:G){
    idx <- sample(1:n, size=j, replace=F)
    Dmean[g] <- calcola_D_norm(XX[idx,],P,0);
  }
  Dmeannj[j,]= quantile(Dmean,prob=c(0.025,0.5,0.975))
}

apply(Dmeannj,2,diff)

#####
Dmeannjold <- matrix(nrow = Nj, ncol = 3) 
Dmeannjold[1,] <- 0


for(j in 2:Nj){
  Dmean <- vector(length=G)
  for(g in 1:G){
    idx <- sample(1:n,size=j,replace=F)
    Dmean[g] <- calcola_D(XX[idx,], P, 0);
  }
  Dmeannjold[j,]=quantile(Dmean,prob=c(0.025,0.5,0.975))
}

par(mfrow=c(1,2))
#distribuzione della distanza dal centroide x campioni di dimensione 1:100
matplot(1:Nj, Dmeannj,type="l")

matplot(1:Nj, Dmeannjold,type="l")

x11()

par(mfrow=c(1,2))
matplot(apply(Dmeannj,2,diff),type="l")
abline(h=0)
matplot(apply(Dmeannjold,2,diff),type="l")
abline(h=0)


Dmeangow <- matrix(nrow=Nj, ncol = 3)
Dmeangow[1,] <- 0

for(j in 2:Nj){
  Dmean <- vector(length=G)
  for(g in 1:G){
    idx <- sample(1:n,size=j,replace=F)
    Dmean[g] <- sum(as.matrix(cluster::daisy(XX[idx,], metric="gower")))
      #calcola_D(XX[idx,], P, 0);
  }
  Dmeangow[j,]=quantile(Dmean,prob=c(0.025,0.5,0.975))
  print(j)
}

par(mfrow=c(1,3))
#distribuzione della distanza dal centroide x campioni di dimensione 1:100
matplot(1:Nj, Dmeannj,type="l")
matplot(1:Nj, Dmeannjold,type="l")
matplot(1:Nj, Dmeangow,type="l")

x11()

par(mfrow=c(1,3))
matplot(apply(Dmeannj,2,diff),type="l")
abline(h=0)
matplot(apply(Dmeannjold,2,diff),type="l")
matplot(apply(Dmeangow,2,diff),type="l")
abline(h=0)

