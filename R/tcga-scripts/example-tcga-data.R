## Install TCGAbiolinks
#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#}
#BiocManager::install("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")

rm(list = ls())
# In the following I am going to follow this tutorial to download data from TCGA:
# https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html
# some lines are slightly different, but the data and procedures are consistent

# Load packages
library("TCGAbiolinks")
library("SummarizedExperiment")

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

lihc_res = getResults(query_TCGA) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table
# head(lihc_res$tissue.definition)
head(lihc_res$sample_type)

summary(as.factor(lihc_res$sample_type))

#we will ignore the small class of recurrent solid tumors. Therefore, we will redo the query as

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

# Download the files from the query. 
# Check that you set your current working directory correctly
#getwd()
GDCdownload(query = query_TCGA)

# Load the actual RNASeq data into R.
tcga_data = GDCprepare(query_TCGA)

dim(tcga_data)

#https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#anatomy-of-a-summarizedexperiment
colnames(colData(tcga_data))

#The table() function (in this context) produces a small summary with the sum of each of the factors present in a given column.
table(tcga_data@colData$vital_status)
table(tcga_data@colData$prior_treatment)
#table(tcga_data@colData$treatment_type)
table(tcga_data@colData$tumor_grade) #not reported...
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)

dim(assay(tcga_data))     # gene expression matrices.
head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples

head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

## Save the data as a file, if you need it later, you can just load this file
## instead of having to run the whole pipeline again
#saveRDS(object = tcga_data,
#        file = "data/tcga_data.RDS",
#        compress = FALSE)
#
#tcga_data <- readRDS("~/Dropbox/PHD/study-treatppmx/data/tcga_data.RDS")