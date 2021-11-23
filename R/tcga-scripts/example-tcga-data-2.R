library(TCGAretriever)

my_csid <- "lgg_ucsf_2014"

blca_pro <- get_genetic_profiles(csid = my_csid)
blca_pro[, 1:2]
