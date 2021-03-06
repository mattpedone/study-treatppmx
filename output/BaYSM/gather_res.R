rm(list=ls())
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_cal_NN.RData")
auxcalnn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_cal_NNIG.RData")
auxcalnnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_coa1_NN.RData")
auxcoa1nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_coa1_NNIG.RData")
auxcoa1nnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_coa2_NN.RData")
auxcoa2nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_aux_coa2_NNIG.RData")
auxcoa2nnig <- t(resPPMX)

load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_cal_NN.RData")
ddcalnn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_cal_NNIG.RData")
ddcalnnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_coa1_NN.RData")
ddcoa1nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_coa1_NNIG.RData")
ddcoa1nnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_coa2_NN.RData")
ddcoa2nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_nu_dd_coa2_NNIG.RData")
ddcoa2nnig <- t(resPPMX)

noupdate <- rbind(auxcalnn, auxcalnnig, auxcoa1nn, auxcoa1nnig, auxcoa2nn, auxcoa2nnig, 
      ddcalnn, ddcalnnig, ddcoa1nn, ddcoa1nnig, ddcoa2nn, ddcoa2nnig)

#save(noupdate, "output/baysm_scenario2/restabnu.RData")

load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_cal_NN.RData")
auxcalnn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_cal_NNIG.RData")
auxcalnnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_coa1_NN.RData")
auxcoa1nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_coa1_NNIG.RData")
auxcoa1nnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_coa2_NN.RData")
auxcoa2nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_aux_coa2_NNIG.RData")
auxcoa2nnig <- t(resPPMX)

load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_cal_NN.RData")
ddcalnn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_cal_NNIG.RData")
ddcalnnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_coa1_NN.RData")
ddcoa1nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_coa1_NNIG.RData")
ddcoa1nnig <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_coa2_NN.RData")
ddcoa2nn <- t(resPPMX)
load("~/Dropbox/PHD/study-treatppmx/output/baysm_scenario2/res_u_dd_coa2_NNIG.RData")
ddcoa2nnig <- t(resPPMX)

update <- rbind(auxcalnn, auxcalnnig, auxcoa1nn, auxcoa1nnig, auxcoa2nn, auxcoa2nnig, 
                  ddcalnn, ddcalnnig, ddcoa1nn, ddcoa1nnig, ddcoa2nn, ddcoa2nnig)

#save(update, "output/baysm_scenario2/restabu.RData")

tab <- cbind(noupdate, update)
rownames(tab) <- c("auxcalnn", " ", "auxcalnnig", " ", "auxcoa1nn", " ", 
"auxcoa1nnig", " ", "auxcoa2nn", " ", "auxcoa2nnig", " ", "ddcalnn", " ", 
"ddcalnnig", " ", "ddcoa1nn", " ", "ddcoa1nnig", " ", "ddcoa2nn", " ", "ddcoa2nnig", " ")

#xtable::xtable(tab, digits = 4)

xtable(tab[-c(1:4, 13:16),c(6:10)])

#xtable::xtable(tab[-c(1:4, 13:16),c(6:10)])
