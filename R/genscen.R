rm(list=ls())
library(treatppmx)
set.seed(121)
### ----
# This script is used to generate all the scenarios
###

# Scenario 1
treatppmx::genmech(npred = 10, progscen = 2, predscen = 1, nset = 30, save = T, 
                   filename = "scenario1")

# Scenario 2
treatppmx::genmech(npred = 25, progscen = 2, predscen = 1, nset = 30, save = T, 
                   filename = "scenario2")

# Scenario 3
treatppmx::genmech(npred = 50, progscen = 2, predscen = 1, nset = 30, save = T, 
                   filename = "scenario3")

# Scenario 4a
treatppmx::genmech(npred = 5, progscen = 2, predscen = 2, nnoise = 5, 
                   nset = 30, save = T, filename = "scenario4a")

# Scenario 4
treatppmx::genmech(npred = 10, progscen = 2, predscen = 2, nnoise = 15, 
                   nset = 30, save = T, filename = "scenario4")

# Scenario 5
treatppmx::genmech(npred = 10, progscen = 2, predscen = 2, nnoise = 40, 
                   nset = 30, save = T, filename = "scenario5")

# Scenario 6
treatppmx::genmech(npred = 20, progscen = 2, predscen = 2, nnoise = 80, 
                   nset = 30, save = T, filename = "scenario6")

# Scenario 1-alt
scenalt1 <- treatppmx::genmech_alt(npred = 10, nset = 30, overlap = 1.0)
save(scenalt1, file = "data/scenalt1.RData")

# Scenario 2-alt
scenalt2 <- treatppmx::genmech_alt(npred = 25, nset = 30, overlap = 1.0)
save(scenalt2, file = "data/scenalt2.RData")

# Scenario 3-alt
scenalt3 <- treatppmx::genmech_alt(npred = 50, nset = 30, overlap = 1.0)
save(scenalt3, file = "data/scenalt3.RData")

# Scenario 4-alt
scenalt4 <- treatppmx::genmech_alt(npred = 10, nset = 30, overlap = .9)
save(scenalt4, file = "data/scenalt4.RData")

# Scenario 5-alt
scenalt5 <- treatppmx::genmech_alt(npred = 25, nset = 30, overlap = .9)
save(scenalt5, file = "data/scenalt5.RData")

# Scenario 6-alt
scenalt6 <- treatppmx::genmech_alt(npred = 50, nset = 30, overlap = .9)
save(scenalt6, file = "data/scenalt6.RData")

# Scenario 7-alt
scenalt7 <- treatppmx::genmech_alt(npred = 10, nset = 30, overlap = .8)
save(scenalt7, file = "data/scenalt7.RData")

# Scenario 8-alt
scenalt8 <- treatppmx::genmech_alt(npred = 25, nset = 30, overlap = .8)
save(scenalt8, file = "data/scenalt8.RData")

# Scenario 9-alt
scenalt9 <- treatppmx::genmech_alt(npred = 50, nset = 30, overlap = .8)
save(scenalt9, file = "data/scenalt9.RData")



