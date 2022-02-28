rm(list=ls())
library(ggplot2)
#library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)

#### ---- Scenarios 1-3 ---- ####
mypath <- c("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-25/leave-one-out/")
load(paste0(mypath, "scen1/mot_ma_pam.rda"))
pam <- HCppCT
load(paste0(mypath, "scen1/mot_ma_km.rda"))
km <- HCppCT
load(paste0(mypath, "scen1/mot_ma_hc.rda"))
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-1/mot.RData")
ppmx <- PPMXCT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mot = value) 

mot1 <- ggplot(df, aes(x=method, y=mot, fill=method)) + 
  ylim(0, 50) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 1") + 
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen2/mot_ma_pam.rda"))
pam <- HCppCT
load(paste0(mypath, "scen2/mot_ma_km.rda"))
km <- HCppCT
load(paste0(mypath, "scen2/mot_ma_hc.rda"))
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-2/mot.RData")
ppmx <- PPMXCT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mot = value) 

mot2 <- ggplot(df, aes(x=method, y=mot, fill=method)) + 
  ylim(0, 50) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 2a") + 
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen3/mot_ma_pam.rda"))
pam <- HCppCT
load(paste0(mypath, "scen3/mot_ma_km.rda"))
km <- HCppCT
load(paste0(mypath, "scen3/mot_ma_hc.rda"))
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-3/mot.RData")
ppmx <- PPMXCT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mot = value) 

mot3 <- ggplot(df, aes(x=method, y=mot, fill=method)) + 
  ylim(0, 50) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 2b") + 
  scale_fill_brewer(palette="Accent") 


load(paste0(mypath, "scen1/mtug_ma_pam.rda"))
pam <- mtug
load(paste0(mypath, "scen1/mtug_ma_km.rda"))
km <- mtug
load(paste0(mypath, "scen1/mtug_ma_hc.rda"))
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-1/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtu = value) 

mtug1 <- ggplot(df, aes(x=method, y=mtu, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen2/mtug_ma_pam.rda"))
pam <- mtug
load(paste0(mypath, "scen2/mtug_ma_km.rda"))
km <- mtug
load(paste0(mypath, "scen2/mtug_ma_hc.rda"))
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-2/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtu = value) 

mtug2 <- ggplot(df, aes(x=method, y=mtu, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen3/mtug_ma_pam.rda"))
pam <- mtug
load(paste0(mypath, "scen3/mtug_ma_km.rda"))
km <- mtug
load(paste0(mypath, "scen3/mtug_ma_hc.rda"))
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-3/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtu = value) 

mtug3 <- ggplot(df, aes(x=method, y=mtu, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")

load(paste0(mypath, "scen1/npc_ma_pam.rda"))
load(paste0(mypath, "scen1/npc_ma_km.rda"))
load(paste0(mypath, "scen1/npc_ma_hc.rda"))
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-1/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc1 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(50, 110) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen2/npc_ma_pam.rda"))
pam <- HCppCUT
load(paste0(mypath, "scen2/npc_ma_km.rda"))
km <- HCppCUT
load(paste0(mypath, "scen2/npc_ma_hc.rda"))
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-2/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc2 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(50, 110) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen3/npc_ma_pam.rda"))
pam <- HCppCUT
load(paste0(mypath, "scen3/npc_ma_km.rda"))
km <- HCppCUT
load(paste0(mypath, "scen3/npc_ma_hc.rda"))
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-3/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc3 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(50, 110) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

p <- gridExtra::grid.arrange(mot1, mot2, mot3, mtug1, mtug2, mtug3, npc1, npc2, npc3, nrow = 3)

ggsave(p, filename = "figs/bp_sim_scen_13.png",  bg = "transparent")

mypath <- c("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-25/leave-one-out/")
load(paste0(mypath, "scen4/mot_ma_pam.rda"))
pam <- HCppCT
load(paste0(mypath, "scen4/mot_ma_km.rda"))
km <- HCppCT
load(paste0(mypath, "scen4/mot_ma_hc.rda"))
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-4/mot.RData")
ppmx <- PPMXCT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mot = value) 

mot1 <- ggplot(df, aes(x=method, y=mot, fill=method)) + 
  ylim(0, 50) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 3a") + 
  scale_fill_brewer(palette="Accent") 

mypath <- c("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-25/leave-one-out/")
load(paste0(mypath, "scen5/mot_ma_pam.rda"))
pam <- HCppCT
load(paste0(mypath, "scen5/mot_ma_km.rda"))
km <- HCppCT
load(paste0(mypath, "scen5/mot_ma_hc.rda"))
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-5/mot.RData")
ppmx <- PPMXCT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mot = value) 

mot2 <- ggplot(df, aes(x=method, y=mot, fill=method)) + 
  ylim(0, 50) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 3b") + 
  scale_fill_brewer(palette="Accent") 


load(paste0(mypath, "scen4/mtug_ma_pam.rda"))
pam <- mtug
load(paste0(mypath, "scen4/mtug_ma_km.rda"))
km <- mtug
load(paste0(mypath, "scen4/mtug_ma_hc.rda"))
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-4/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtu = value) 

mtug1 <- ggplot(df, aes(x=method, y=mtu, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen5/mtug_ma_pam.rda"))
pam <- mtug
load(paste0(mypath, "scen5/mtug_ma_km.rda"))
km <- mtug
load(paste0(mypath, "scen5/mtug_ma_hc.rda"))
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-5/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtu = value) 

mtug2 <- ggplot(df, aes(x=method, y=mtu, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen4/npc_ma_pam.rda"))
pam <- HCppCUT
load(paste0(mypath, "scen4/npc_ma_km.rda"))
km <- HCppCUT
load(paste0(mypath, "scen4/npc_ma_hc.rda"))
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-4/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc1 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(50, 110) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(paste0(mypath, "scen5/npc_ma_pam.rda"))
pam <- HCppCUT
load(paste0(mypath, "scen5/npc_ma_km.rda"))
km <- HCppCUT
load(paste0(mypath, "scen5/npc_ma_hc.rda"))
hc <- HCPPCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/sigma-01/leave-one-out/scen-5/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc2 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(50, 110) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

p <- gridExtra::grid.arrange(mot1, mot2, mtug1, mtug2, npc1, npc2, nrow = 3)

ggsave(p, filename = "figs/bp_sim_scen_45.png",  bg = "transparent")

