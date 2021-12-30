rm(list=ls())
##Scenario 1
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_hc.RData")
mot <- c(mean(MOT), sd(MOT))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_hc.RData")
mtug <- c(mean(MTUg), sd(MTUg))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_hc.RData")
npc <- c(mean(NPC), sd(NPC))
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_km.RData")
mot <- c(mean(MOT), sd(MOT))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_km.RData")
mtug <- c(mean(MTUg), sd(MTUg))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_km.RData")
npc <- c(mean(NPC), sd(NPC))
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_pam.RData")
mot <- c(mean(MOT), sd(MOT))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_pam.RData")
mtug <- c(mean(MTUg), sd(MTUg))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_pam.RData")
npc <- c(mean(NPC), sd(NPC))
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/res.RData")
ppmx <- t(resPPMX[1:3,])

rbind(pam, km, hc, ppmx)

##Scenario 2
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_hc.RData")
mot <- c(mean(MOT), sd(MOT))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_hc.RData")
mtug <- c(mean(MTUg), sd(MTUg))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_hc.RData")
npc <- c(mean(NPC), sd(NPC))
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_km.RData")
mot <- c(mean(MOT), sd(MOT))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_km.RData")
mtug <- c(mean(MTUg), sd(MTUg))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_km.RData")
npc <- c(mean(NPC), sd(NPC))
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

##pam
#load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_pam.RData")
#mot <- c(mean(MOT), sd(MOT))
#load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_pam.RData")
#mtug <- c(mean(MTUg), sd(MTUg))
#load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_pam.RData")
#npc <- c(mean(NPC), sd(NPC))
#pam <- t(rbind(mot, mtug, npc))
#rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/res.RData")
ppmx <- t(resPPMX[1:3,])

rbind(km, hc, ppmx)
