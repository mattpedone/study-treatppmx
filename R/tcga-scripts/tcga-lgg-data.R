rm(list = ls())
# Load packages
library("TCGAbiolinks")
library("SummarizedExperiment")

TCGAbiolinks:::getProjectSummary("TCGA-LGG")

query_lgg = GDCquery(
  project = "TCGA-LGG",
  data.category = "Proteome Profiling",
  #data.type = "miRNA Expression Quantification",
  #sample.type = "Primary Tumor", 
  #experimental.strategy = "miRNA-Seq",
  #workflow.type = "BCGSC miRNA Profiling",
  #platform = "RPPA",
  legacy = FALSE)

lgg_res = getResults(query_lgg) 
colnames(lgg_res) 

# Commento Query
### ---- 
# Non sono sicuro della query fatta. data category è corretto? Sono indeciso tra
# Transcriptome Profiling e Proteome Profiling.
# Se scelgo Proteome Profiling l'unica experimental.strategy è Reverse Phase Protein Array,
# ma non riesco a specificare un workflow.type valido

# Se metto 
#data.category = c("Transcriptome Profiling"),
# ci vuole circa 10 min.
# Sulla base di quanto segue credo non sia corretto
# > table(lgg_res$type)
# 
# gene_expression mirna_expression 
# 1587             1060 
# > table(lgg_res$data_type)
# 
# Gene Expression Quantification Isoform Expression Quantification 
# 1587                               530 
# miRNA Expression Quantification 
# 530 
# > table(lgg_res$experimental_strategy)
# 
# miRNA-Seq   RNA-Seq 
# 1060      1587 
# > table(lgg_res$analysis_workflow_type)
# 
# BCGSC miRNA Profiling        HTSeq - Counts          HTSeq - FPKM 
# 1060                   529                   529 
# HTSeq - FPKM-UQ 
# 529 
### ---- 
# Fine commento query

table(lgg_res$sample_type)
# Maybe we need to ignore the class of Recurrent Tumors. The query would be:
#query_TCGA = GDCquery(
#  project = "TCGA-LGG",
#  data.category = c("Proteome Profiling"),
#  sample.type = c("Primary Tumor"))

# Download the files from the query. 
# Check that you set your current working directory correctly
getwd()
GDCdownload(query = query_lgg)

# Load the data into R.
query_lgg
lgg_data = GDCprepare(query_lgg, summarizedExperiment = T)

dim(lgg_data)

#https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#anatomy-of-a-summarizedexperiment
colnames(colData(lgg_data))

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
