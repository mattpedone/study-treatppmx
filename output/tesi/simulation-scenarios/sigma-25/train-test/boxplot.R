rm(list=ls())
library(ggplot2)
#library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)

load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_hc.RData")
hc <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_km.RData")
km <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot_pam.RData")
pam <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mot.RData")
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

load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_hc.RData")
hc <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot_km.RData")
km <- MOT
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mot.RData")
ppmx <- PPMXCT

df <- tibble(km = km, hc = hc, ppmx = ppmx)

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


load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_hc.RData")
hc <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_km.RData")
km <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug_pam.RData")
pam <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/mtug.RData")
ppmx <- PPMXpp

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug1 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 1") + 
  scale_fill_brewer(palette="Accent") 

load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_hc.RData")
hc <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug_km.RData")
km <- MTUg
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/mtug.RData")
ppmx <- PPMXpp

df <- tibble(km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, mtug = value) 

mtug2 <- ggplot(df, aes(x=method, y=mtug, fill=method)) + 
  ylim(0, 1) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 2") + 
  scale_fill_brewer(palette="Accent") 

load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_hc.RData")
hc <- NPC
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_km.RData")
km <- NPC
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc_pam.RData")
pam <- NPC
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-1/npc.RData")
ppmx <- PPMXCUT

df <- tibble(pam = pam, km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc1 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(0, 20) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 1") + 
  scale_fill_brewer(palette="Accent") 

load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_hc.RData")
hc <- NPC
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc_km.RData")
km <- NPC
load("~/Dropbox/PHD/study-treatppmx/output/simulation-scenarios/train-test/scen-alt-2/npc.RData")
ppmx <- PPMXCUT

df <- tibble(km = km, hc = hc, ppmx = ppmx)

df <- df %>%
  melt  %>%
  rename(method = variable, npc = value) 

npc2 <- ggplot(df, aes(x=method, y=npc, fill=method)) + 
  ylim(0, 20) +
  geom_boxplot(color="black", fill="grey", alpha=0.3) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle("Scenario 2") + 
  scale_fill_brewer(palette="Accent") 


gridExtra::grid.arrange(mot1, mot2, mtug1, mtug2, npc1, npc2, nrow = 3)

#ggsave(p, filename = "figs/bp_sim_scen_46.png",  bg = "transparent")
