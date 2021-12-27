rm(list=ls())
load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/res_ma_pam.rda"))
pam <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/res_ma_km.rda")
km <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/res_ma_hc.rda")
hc <- resHCpp
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/res.RData")
ppmx <- resPPMX

scen6 <- rbind(t(pam), t(km), t(hc), t(ppmx[1:3,]))
rownames(scen6) <- c("pam", "", "km", "", "hc", "", "ppmx", "")

save(scen6, file = "output/simulation-scenarios/results/scen6.RData")
