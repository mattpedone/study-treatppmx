rm(list=ls())
library(ggplot2)
#library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)

#### ---- Scenarios 1-3 ---- ####

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mot.RData")
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

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mot.RData")
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
  ggtitle("Scenario 2") + 
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mot.RData")
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
  ggtitle("Scenario 3") + 
  scale_fill_brewer(palette="Accent") 


load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug1 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug2 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug3 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")


load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen1/npc.RData")
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

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen2/npc.RData")
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

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen3/npc.RData")
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

#### ---- Scenarios 4-6 ---- ####

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mot.RData")
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
  ggtitle("Scenario 4") + 
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mot.RData")
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
  ggtitle("Scenario 5") + 
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mot_ma_pam.rda"))
pam <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mot_ma_km.rda")
km <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mot_ma_hc.rda")
hc <- HCppCT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mot.RData")
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
  ggtitle("Scenario 6") + 
  scale_fill_brewer(palette="Accent") 


load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug1 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug2 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") 

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mtug_ma_pam.rda"))
pam <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mtug_ma_km.rda")
km <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mtug_ma_hc.rda")
hc <- mtug
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug3 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0.25, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")


load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen4/npc.RData")
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

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen5/npc.RData")
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

load(("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/npc_ma_pam.rda"))
pam <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/npc_ma_km.rda")
km <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/npc_ma_hc.rda")
hc <- HCppCUT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/scen6/npc.RData")
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

ggsave(p, filename = "figs/bp_sim_scen_46.png",  bg = "transparent")
