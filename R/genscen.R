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
scenalt1 <- treatppmx::genmech_alt(npred = 25, nset = 30, overlap = .8)
save(scenalt1, file = "data/scenalt1.RData")

# Scenario 2-alt
scenalt2 <- treatppmx::genmech_alt(npred = 50, nset = 30, overlap = .8)
save(scenalt2, file = "data/scenalt2.RData")

# Scenario 3-alt
scenalt3 <- treatppmx::genmech_alt(npred = 25, nset = 30, overlap = .9)
save(scenalt3, file = "data/scenalt3.RData")

# Scenario 4-alt
scenalt4 <- treatppmx::genmech_alt(npred = 50, nset = 30, overlap = .9)
save(scenalt4, file = "data/scenalt4.RData")

