rm(list=ls())
##Scenario 1
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt1 <- rbind(pam, km, hc, ppmx)

##Scenario 2
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt2 <- rbind(pam, km, hc, ppmx)

##Scenario 3
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-3/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt3 <- rbind(pam, km, hc, ppmx)

##Scenario 4
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-4/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt4 <- rbind(pam, km, hc, ppmx)

##Scenario 5
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-5/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt5 <- rbind(pam, km, hc, ppmx)

##Scenario 6
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-6/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt6 <- rbind(pam, km, hc, ppmx)

##Scenario 7
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-7/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt7 <- rbind(pam, km, hc, ppmx)

##Scenario 8
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-8/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt8 <- rbind(pam, km, hc, ppmx)

##Scenario 9
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-9/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt9 <- rbind(pam, km, hc, ppmx)

##Scenario 10
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-10/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt10 <- rbind(pam, km, hc, ppmx)

##Scenario 11
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-11/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt11 <- rbind(pam, km, hc, ppmx)

##Scenario 12
#hc
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mot_hc.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mtug_hc.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/npc_hc.RData")
npc <- NPC
hc <- t(rbind(mot, mtug, npc))
rownames(hc) <- c("mean", "sd")

#km
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mot_km.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mtug_km.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/npc_km.RData")
npc <- NPC
km <- t(rbind(mot, mtug, npc))
rownames(km) <- c("mean", "sd")

#pam
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mot_pam.RData")
mot <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/mtug_pam.RData")
mtug <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/npc_pam.RData")
npc <- NPC
pam <- t(rbind(mot, mtug, npc))
rownames(pam) <- c("mean", "sd")

#ppmx
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-12/res.RData")
ppmx <- t(resPPMX[1:3,])

scenalt12 <- rbind(pam, km, hc, ppmx)