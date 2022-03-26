rm(list=ls())

library(dplyr)
library(tibble)
library(tidyselect)
library(ggplot2)
library(forcats)
library(ggpubr)

dir <- "~/Dropbox/PHD/study-treatppmx/output/journal/simulations"

#MOT A
load(paste0(dir, "/ma/scen-alt-2/mot_hc.RData"))
hc1 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-2/mot_km.RData"))
km1 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-2/mot_pam.RData"))
pam1 <- MOT[1]
load(paste0(dir, "/scen2mot.RData"))
ppmx1 <- mean(PPMXCT)

load(paste0(dir, "/ma/scen-alt-5/mot_hc.RData"))
hc2 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-5/mot_km.RData"))
km2 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-5/mot_pam.RData"))
pam2 <- MOT[1]
load(paste0(dir, "/scen5mot.RData"))
ppmx2 <- mean(PPMXCT)


load(paste0(dir, "/ma/scen-alt-8/mot_hc.RData"))
hc3 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-8/mot_km.RData"))
km3 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-8/mot_pam.RData"))
pam3 <- MOT[1]
load(paste0(dir, "/scen8mot.RData"))
ppmx3 <- mean(PPMXCT)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(mot=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(mot=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(mot=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(mot=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("a", nrow(res)))
mota <- ggplot(res, aes(x=heterogeneity, y=mot, group=method, color=method)) + 
  ylim(2, 9) +
  ylab("MOT") +
  xlab("") +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") + 
  facet_grid(.~scen)

#MOT B
load(paste0(dir, "/ma/scen-alt-3/mot_hc.RData"))
hc1 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-3/mot_km.RData"))
km1 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-3/mot_pam.RData"))
pam1 <- MOT[1]
load(paste0(dir, "/scen3mot.RData"))
ppmx1 <- mean(PPMXCT)

load(paste0(dir, "/ma/scen-alt-6/mot_hc.RData"))
hc2 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-6/mot_km.RData"))
km2 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-6/mot_pam.RData"))
pam2 <- MOT[1]
load(paste0(dir, "/scen6mot.RData"))
ppmx2 <- mean(PPMXCT)


load(paste0(dir, "/ma/scen-alt-9/mot_hc.RData"))
hc3 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-9/mot_km.RData"))
km3 <- MOT[1]
load(paste0(dir, "/ma/scen-alt-9/mot_pam.RData"))
pam3 <- MOT[1]
load(paste0(dir, "/scen9mot.RData"))
ppmx3 <- mean(PPMXCT)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(mot=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(mot=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(mot=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(mot=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("b", nrow(res)))
motb <- ggplot(res, aes(x=heterogeneity, y=mot, group=method, color=method)) +
  ylim(2, 9) +
  ylab("MOT") + 
  xlab("") +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")  + 
  facet_grid(.~scen)

#MTU A
load(paste0(dir, "/ma/scen-alt-2/mtug_hc.RData"))
hc1 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-2/mtug_km.RData"))
km1 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-2/mtug_pam.RData"))
pam1 <- MTUg[1]
load(paste0(dir, "/scen2mtug.RData"))
ppmx1 <- mean(PPMXpp)

load(paste0(dir, "/ma/scen-alt-5/mtug_hc.RData"))
hc2 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-5/mtug_km.RData"))
km2 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-5/mtug_pam.RData"))
pam2 <- MTUg[1]
load(paste0(dir, "/scen5mtug.RData"))
ppmx2 <- mean(PPMXpp)


load(paste0(dir, "/ma/scen-alt-8/mtug_hc.RData"))
hc3 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-8/mtug_km.RData"))
km3 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-8/mtug_pam.RData"))
pam3 <- MTUg[1]
load(paste0(dir, "/scen8mtug.RData"))
ppmx3 <- mean(PPMXpp)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(mtu=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(mtu=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(mtu=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(mtu=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("a", nrow(res)))
mtua <- ggplot(res, aes(x=heterogeneity, y=mtu, group=method, color=method)) +
  ylim(0.45, 0.8) +
  labs( x = "", y = expression(paste("%", Delta , MTU[l]))) +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")  + 
  facet_grid(.~scen)

#MOT B
load(paste0(dir, "/ma/scen-alt-3/mtug_hc.RData"))
hc1 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-3/mtug_km.RData"))
km1 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-3/mtug_pam.RData"))
pam1 <- MTUg[1]
load(paste0(dir, "/scen3mtug.RData"))
ppmx1 <- mean(PPMXpp)

load(paste0(dir, "/ma/scen-alt-6/mtug_hc.RData"))
hc2 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-6/mtug_km.RData"))
km2 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-6/mtug_pam.RData"))
pam2 <- MTUg[1]
load(paste0(dir, "/scen6mtug.RData"))
ppmx2 <- mean(PPMXpp)


load(paste0(dir, "/ma/scen-alt-9/mtug_hc.RData"))
hc3 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-9/mtug_km.RData"))
km3 <- MTUg[1]
load(paste0(dir, "/ma/scen-alt-9/mtug_pam.RData"))
pam3 <- MTUg[1]
load(paste0(dir, "/scen9mtug.RData"))
ppmx3 <- mean(PPMXpp)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(mtu=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(mtu=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(mtu=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(mtu=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("b", nrow(res)))
mtub <- ggplot(res, aes(x=heterogeneity, y=mtu, group=method, color=method)) +
  ylim(0.45, 0.8) +
  labs( x = "", y = expression(paste("% ", Delta, MTU[l]))) +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")  + 
  facet_grid(.~scen)

#NPC A
load(paste0(dir, "/ma/scen-alt-2/npc_hc.RData"))
hc1 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-2/npc_km.RData"))
km1 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-2/npc_pam.RData"))
pam1 <- NPC[1]
load(paste0(dir, "/scen2npc.RData"))
ppmx1 <- mean(PPMXCUT)

load(paste0(dir, "/ma/scen-alt-5/npc_hc.RData"))
hc2 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-5/npc_km.RData"))
km2 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-5/npc_pam.RData"))
pam2 <- NPC[1]
load(paste0(dir, "/scen5npc.RData"))
ppmx2 <- mean(PPMXCUT)


load(paste0(dir, "/ma/scen-alt-8/npc_hc.RData"))
hc3 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-8/npc_km.RData"))
km3 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-8/npc_pam.RData"))
pam3 <- NPC[1]
load(paste0(dir, "/scen8npc.RData"))
ppmx3 <- mean(PPMXCUT)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(npc=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(npc=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(npc=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(npc=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("a", nrow(res)))
npca <- ggplot(res, aes(x=heterogeneity, y=npc, group=method, color=method)) +
  ylim(5, 18) +
  xlab("") +
  ylab("NPC") +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent")  + 
  facet_grid(.~scen)

#NPC B
load(paste0(dir, "/ma/scen-alt-3/npc_hc.RData"))
hc1 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-3/npc_km.RData"))
km1 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-3/npc_pam.RData"))
pam1 <- NPC[1]
load(paste0(dir, "/scen3npc.RData"))
ppmx1 <- mean(PPMXCUT)

load(paste0(dir, "/ma/scen-alt-6/npc_hc.RData"))
hc2 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-6/npc_km.RData"))
km2 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-6/npc_pam.RData"))
pam2 <- NPC[1]
load(paste0(dir, "/scen6npc.RData"))
ppmx2 <- mean(PPMXCUT)


load(paste0(dir, "/ma/scen-alt-9/npc_hc.RData"))
hc3 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-9/npc_km.RData"))
km3 <- NPC[1]
load(paste0(dir, "/ma/scen-alt-9/npc_pam.RData"))
pam3 <- NPC[1]
load(paste0(dir, "/scen9npc.RData"))
ppmx3 <- mean(PPMXCUT)

hf <- as.factor(c("1", "2", "3"))
hc <- data.frame(npc=rbind(hc1, hc2, hc3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("hc", 3))
km <- data.frame(npc=rbind(km1, km2, km3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("km", 3))
pam <- data.frame(npc=rbind(pam1, pam2, pam3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("pam", 3))
ppmx <- data.frame(npc=rbind(ppmx1, ppmx2, ppmx3), heterogeneity = fct_relevel(hf, "1", "2", "3"), method = rep("ppmx", 3))

res <- rbind(hc, km, pam, ppmx)
res <- cbind(res, scen=rep("b", nrow(res)))
npcb <- ggplot(res, aes(x=heterogeneity, y=npc, group=method, color=method)) +
  ylim(5, 18) +
  xlab("") +
  ylab("NPC") +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Accent") + 
  facet_grid(.~scen)

#p <- ggarrange(mota, mtua, npca, motb, mtub, npcb, nrow=2, ncol = 3, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
p <- ggarrange(mota, motb, mtua, mtub, npca, npcb, nrow=3, ncol = 2, common.legend = TRUE, legend="bottom")#, panel.border = element_blank())
p
#ggsave(p, device = "pdf", path = "figs", filename = "trend_het.pdf")
