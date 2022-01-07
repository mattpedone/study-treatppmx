#### ---- table \kappa sigma LGG ---- ####
rm(list=ls())
load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/res_kappasigma.RData")
kappasigma <- res

ks_nn <- rbind(cbind(kappasigma[[1]], kappasigma[[2]], kappasigma[[3]]), 
               cbind(kappasigma[[4]], kappasigma[[5]], kappasigma[[6]]), 
               cbind(kappasigma[[7]], kappasigma[[8]], kappasigma[[9]]))

#colnames(ks_nn) <- c("0.01", "", "0.05", "", "0.25", "")
#rownames(ks_nn) <- c("1", "", "", "", "", "10", "", "", "", "", "40", "", "", "", "")

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/res_kappasigma.RData")
kappasigma <- res

ks_nnig <- rbind(cbind(kappasigma[[1]], kappasigma[[2]], kappasigma[[3]]), 
                 cbind(kappasigma[[4]], kappasigma[[5]], kappasigma[[6]]), 
                 cbind(kappasigma[[7]], kappasigma[[8]], kappasigma[[9]]))

#colnames(ks_nnig) <- c("0.01", "", "0.05", "", "0.25", "")
#rownames(ks_nnig) <- c("1", "", "", "", "", "10", "", "", "", "", "40", "", "", "", "")

ks <- cbind(ks_nn, ks_nnig)
# results in terms of MOT and esm
# kappa values on the rows, sigma value on the columns
#first 3 columns are DDNN, latter 3 DDNNIG
ks

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/clu_kappasigma.RData")
kappasigma <- clu

ks_nn <- rbind(cbind(kappasigma[[1]], kappasigma[[2]], kappasigma[[3]]), 
               cbind(kappasigma[[4]], kappasigma[[5]], kappasigma[[6]]), 
               cbind(kappasigma[[7]], kappasigma[[8]], kappasigma[[9]]))
ks_nn <- ks_nn[-c(2, 5, 8),]

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/clu_kappasigma.RData")
kappasigma <- clu

ks_nnig <- rbind(cbind(kappasigma[[1]], kappasigma[[2]], kappasigma[[3]]), 
                 cbind(kappasigma[[4]], kappasigma[[5]], kappasigma[[6]]), 
                 cbind(kappasigma[[7]], kappasigma[[8]], kappasigma[[9]]))
ks_nnig <- ks_nnig[-c(2, 5, 8),]
ks_clu <- cbind(ks_nn, ks_nnig)

ks_clu

#### ---- table sigma S0 LGG ---- ####
rm(list=ls())
load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/res_sigmas0.RData")
sigmas0 <- res

ss_nn <- rbind(cbind(sigmas0[[1]], sigmas0[[2]], sigmas0[[3]]), 
               cbind(sigmas0[[4]], sigmas0[[5]], sigmas0[[6]]), 
               cbind(sigmas0[[7]], sigmas0[[8]], sigmas0[[9]]))

colnames(ss_nn) <- c(".1", "", "1", "", "10", "")
rownames(ss_nn) <- c("1", "", "", "", "", "10", "", "", "", "", "50", "", "", "", "")

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/res_sigmas0.RData")
sigmas0 <- res

ss_nnig <- rbind(cbind(sigmas0[[1]], sigmas0[[2]], sigmas0[[3]]), 
                 cbind(sigmas0[[4]], sigmas0[[5]], sigmas0[[6]]), 
                 cbind(sigmas0[[7]], sigmas0[[8]], sigmas0[[9]]))

colnames(ss_nnig) <- c(".1", "", "1", "", "10", "")
rownames(ss_nnig) <- c("1", "", "", "", "", "10", "", "", "", "", "50", "", "", "", "")

ss <- cbind(ss_nn, ss_nnig)

ss

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/clu_sigmas0.RData")
sigmas0 <- clu

ss_nn <- rbind(cbind(sigmas0[[1]], sigmas0[[2]], sigmas0[[3]]), 
               cbind(sigmas0[[4]], sigmas0[[5]], sigmas0[[6]]), 
               cbind(sigmas0[[7]], sigmas0[[8]], sigmas0[[9]]))

ss_nn <- ss_nn[-c(2, 5, 8),]

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/clu_sigmas0.RData")
sigmas0 <- clu

ss_nnig <- rbind(cbind(sigmas0[[1]], sigmas0[[2]], sigmas0[[3]]), 
                 cbind(sigmas0[[4]], sigmas0[[5]], sigmas0[[6]]), 
                 cbind(sigmas0[[7]], sigmas0[[8]], sigmas0[[9]]))

ss_nnig <- ss_nnig[-c(2, 5, 8),]

ss <- cbind(ss_nn, ss_nnig)

ss

#### ---- table sigma20 LGG ---- ####
rm(list=ls())
load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/res_sigma20.RData")
sigma <- res

s_nn <- rbind(cbind(sigma[[1]], sigma[[2]], sigma[[3]]))

colnames(s_nn) <- c("1", "", "2", "", "10", "")

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/res_sigma20.RData")
sigma <- res

s_nnig <- rbind(cbind(sigma[[1]], sigma[[2]], sigma[[3]]))

colnames(s_nnig) <- c("1", "", "2", "", "10", "")

s <- cbind(s_nn, s_nnig)

s

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNN/scen1/clu_sigma20.RData")
sigma <- clu

s_nn <- rbind(cbind(sigma[[1]], sigma[[2]], sigma[[3]]))
s_nn <- s_nn[-c(2),]

load("~/Dropbox/PHD/study-treatppmx/output/simulation-sensitivity/DDNNIG/scen1/clu_sigma20.RData")
sigma <- clu

s_nnig <- rbind(cbind(sigma[[1]], sigma[[2]], sigma[[3]]))
s_nnig <- s_nnig[-c(2),]

s_clu <- cbind(s_nn, s_nnig)
