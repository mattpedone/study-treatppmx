rm(list=ls())
load("output/simulation-scenarios/results/scen1.RData")
load("output/simulation-scenarios/results/scen2.RData")
load("output/simulation-scenarios/results/scen3.RData")
load("output/simulation-scenarios/results/scen4.RData")
load("output/simulation-scenarios/results/scen5.RData")


xtable::xtable(cbind(scen1, scen2, scen3), digits = 4)

xtable::xtable(cbind(scen4, scen5), digits = 4)
