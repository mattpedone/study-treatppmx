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
scenalt1 <- treatppmx::genmech_alt(npred = 10, nset = 30)
save(scenalt1, file = "data/scenalt1.RData")

# Scenario 2-alt
scenalt2 <- treatppmx::genmech_alt(npred = 50, nset = 30)
save(scenalt2, file = "data/scenalt2.RData")

# Scenario 1-alt-bis
scenalt1b <- treatppmx::genmech_alt(npred = 10, nset = 30)
save(scenalt1b, file = "data/scenalt1b.RData")

# Scenario 2-alt-bis
scenalt2b <- treatppmx::genmech_alt(npred = 50, nset = 30)
save(scenalt2b, file = "data/scenalt2b.RData")
