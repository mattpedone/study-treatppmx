rm(list=ls())
library(treatppmx)

### ----
# This script is used to generate all the scenarios
###

# this is what I called scenario 2, but with 20 covariates instead of 10
treatppmx::genmech(npred = 20, progscen = 1, predscen = 1, nset = 30, save = T, 
                   filename = "modscenario2")
