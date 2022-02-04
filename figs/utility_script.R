rm(list=ls())
library(ggplot2)

load("data/scenario1.rda")
X <- data.frame(mydata)
Z <- data.frame(cbind(myz2, myz3))
k <- 1
Y <- mytot[,,k]
#idx <- sort(sample(1:nrow(Y), 1, replace = F))
#print(idx)
wk <- c(0, 40, 100)
df <- data.frame(sort(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk)))
myorder <- c(sample(1:117, 117, FALSE), sample(118:152, 35, FALSE))
df1 <- data.frame(df[myorder,])
colnames(df1) <- c("Difference in utility")

load("data/scenario2.rda")
X <- data.frame(mydata)
Z <- data.frame(cbind(myz2, myz3))
k <- 1
Y <- mytot[,,k]
#idx <- sort(sample(1:nrow(Y), 1, replace = F))
#print(idx)
wk <- c(0, 40, 100)
df <- data.frame(sort(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk)))
myorder <- c(sample(1:110, 110, FALSE), sample(111:152, 42, FALSE))
df2 <- data.frame(df[myorder,])
colnames(df2) <- c("Difference in utility")

load("data/scenario3.rda")
X <- data.frame(mydata)
Z <- data.frame(cbind(myz2, myz3))
k <- 1
Y <- mytot[,,k]
#idx <- sort(sample(1:nrow(Y), 1, replace = F))
#print(idx)
wk <- c(0, 40, 100)
df <- data.frame(sort(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk)))
myorder <- c(sample(1:109, 109, FALSE), sample(110:152, 43, FALSE))
df3 <- data.frame(df[myorder,])
colnames(df3) <- c("Difference in utility")

df4 <- data.frame(rep(0, 152))
colnames(df4) <- c("Difference in utility")

val = as.matrix(rbind(df1, df2, df3, df4), ncol = 1)

df <- data.frame(x=rep(c(1:152),4), val=val, 
                 variable=rep(c("Scenario 1","Scenario2","Scenario3", "Reference"), 
                              each=152))

#df <- cbind(x = as.numeric(row.names(df)), df)
p <- ggplot2::ggplot(df) + 
  geom_line(aes(x = x, y = val, colour = variable, linetype = variable)) + theme_classic() + #theme(legend.title="") + 
  #geom_line(aes(x = x, y = val, linetype = variable)) + #theme_classic() + #theme(legend.title="") + 
  #geom_line(aes(x = x, y = val, colour = variable)) + theme_classic() + #theme(legend.title="") + 
  xlab("Patients")  +  ylab("Difference in utility") +
  scale_x_continuous(breaks = round(seq(min(df[1]), max(df[1]), by = 10),0)) +
  scale_y_continuous(breaks = round(seq(min(df[2]), max(df[2]), by = 12),1)) +
  #scale_colour_manual(values = c(4, RColorBrewer::brewer.pal(9, "Greys")[c(4, 6, 9)]))
  scale_color_manual(values = c("#F8766D","#E69F00", "#56B4E9", "#009E73"))
p <- p + theme(legend.title = element_blank())
p

ggsave(p, filename = "figs/utility13.png",  bg = "transparent")
