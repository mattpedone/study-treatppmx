rm(list = ls())
# Load packages
## pckgs for proteins & clinical covariates
library("TCGAbiolinks")
library("SummarizedExperiment")

## pckgs for data wrangling
library(dplyr)
library(tidyverse)
library(data.table)
library(tibble)
library(tidyfst)
library(stringr)

# Available data for Lower Grade Glioma
TCGAbiolinks:::getProjectSummary("TCGA-LGG", legacy = T)

# Proteins Tibble
query_protein = GDCquery(
  project = "TCGA-LGG",
  data.category = "Protein expression",
  legacy = T,
  data.format = "BCR Biotab")

GDCdownload(query = query_protein)
lgg_protein = GDCprepare(query_protein, summarizedExperiment = T)

colnames(lgg_protein) <- substring(colnames(lgg_protein), first = 1, last=12)
lgg_protein <- lgg_protein[,as.logical(1-duplicated(colnames(lgg_protein))),]

protein <- data.table::transpose(lgg_protein)
colnames(protein) <- rownames(lgg_protein)
rownames(protein) <- colnames(lgg_protein)

protein <- protein %>%
  janitor::row_to_names(1)

# Clinical Tibbles
query_clinical <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical", file.type = "xml")

GDCdownload(query_clinical)

clinical <- GDCprepare_clinic(query_clinical, clinical.info = "patient")
clinical_drug <- GDCprepare_clinic(query_clinical, clinical.info = "drug")
clinical_radiation <- GDCprepare_clinic(query_clinical, clinical.info = "radiation")

# I need to remove follow-up records. I keep only the first.
# First, I select only the columns I am interested in, then I select the older record

##Clinical
clinical_nfu <- clinical %>%
  select(bcr_patient_barcode, gender, days_to_birth, race_list, 
         has_drugs_information, has_radiations_information)

rep <- clinical_nfu %>% 
  group_by(bcr_patient_barcode) %>% 
  filter(n()>1)

# All the patients that have multiple records are duplicated records collected 
# in the same moment that carry the very same information. I will keep the first one

clinical_nfu <- clinical_nfu %>% 
  distinct()

##Clinical_drug
clinical_drug_nfu <- clinical_drug %>%
  select(bcr_patient_barcode, therapy_types)

# For the patients that have duplicated records I keep the first one

clinical_drug_nfu <- clinical_drug_nfu %>% 
  distinct()

# For those that have received a combination of therapies I set the therapy as "advanced"

rep <- clinical_drug_nfu %>% 
  group_by(bcr_patient_barcode) %>% 
  filter(n()>1)

dup_id <- unique(rep$bcr_patient_barcode)

sing <- clinical_drug_nfu %>% 
  group_by(bcr_patient_barcode) %>% 
  filter(n()<2)

sing_id <- sing$bcr_patient_barcode

treatment <- data.frame(bcr_patient_barcode = c(dup_id, sing_id))

#se la combinazione non prevede Targeted Molecular therapy non è advanced
treat_dup <- data.frame(id = dup_id, treatment = rep("Standard", length(dup_id)))
idtmt <- subset(rep, rep$therapy_types == "Targeted Molecular therapy")$bcr_patient_barcode
for(i in 1:nrow(treat_dup)){
  if(treat_dup$id[1] %in% idtmt){
    treat_dup$treatment <- "Advanced"
  }
}
treatment <- cbind(treatment, treatment = c(treat_dup$treatment, 
                                            as.character(sing$therapy_types)))

clinical_nfu <- left_join(clinical_nfu, treatment)

clinical_nfu <- clinical_nfu %>%
  mutate(treatment=replace(treatment, has_radiations_information == "YES", "Advanced")) %>%
  mutate(treatment=replace(treatment, treatment == "Immunotherapy", "Standard")) %>%
  mutate(treatment=replace(treatment, treatment == "Chemotherapy", "Standard")) %>%
  as.data.frame() 

clinical_nfu <- subset(clinical_nfu, clinical_nfu$treatment != "<NA>")

table(clinical_nfu$treatment)

#outcome
clinical$primary_therapy_outcome_success

clinical_radiation$measure_of_response

clinical_drug$measure_of_response




Reduce(intersect, list(rownames(protein),clinical$bcr_patient_barcode))
       
       , 
                       clinical_drug$bcr_patient_barcode, clinical_radiation$bcr_patient_barcode))
