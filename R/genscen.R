rm(list=ls())
library(treatppmx)
set.seed(121)
### ----
# This script is used to generate all the scenarios
###

# this is scenario 2
treatppmx::genmech(npred = 10, progscen = 1, predscen = 1, nset = 30, save = T, 
                   filename = "scenario2")

# this is what I called modified scenario 2, but with 20 covariates instead of 10
treatppmx::genmech(npred = 20, progscen = 1, predscen = 1, nset = 30, save = T, 
                   filename = "modscenario2")

# this is scenario 3
treatppmx::genmech(npred = 10, progscen = 2, predscen = 1, nset = 30, save = T, 
                   filename = "scenario3")

# this is scenario 4
treatppmx::genmech(npred = 25, progscen = 2, predscen = 1, nset = 30, save = T, 
                   filename = "scenario4")
