rm(list=ls())
load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/res_ma_pam.rda"))
pam <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/res_ma_km.rda")
km <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/res_ma_hc.rda")
hc <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/res.RData")
ppmx <- resPPMX

scen5 <- rbind(t(pam), t(km), t(hc), t(ppmx[1:3,]))
rownames(scen5) <- c("pam", "", "km", "", "hc", "", "ppmx", "")

save(scen5, file = "output/simulation-scenarios/results/scen5.RData")