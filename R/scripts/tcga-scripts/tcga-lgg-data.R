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
TCGAbiolinks:::getProjectSummary("TCGA-LGG")

# Proteins Tibble
query_lgg = GDCquery(
  project = "TCGA-LGG",
  data.category = "Proteome Profiling",
  #sample.type = "Primary Tumor", 
  data.format = "BCR Biotab",
  legacy = FALSE)

GDCdownload(query = query_lgg)
lgg_data = GDCprepare(query_lgg, summarizedExperiment = T)

# Clinical Tibble
query_lgg_clin = GDCquery(
  project = "TCGA-LGG",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
#  sample.type = "Primary Tumor", 
  data.format = "BCR Biotab")

GDCdownload(query = query_lgg_clin)
lgg_data_clin = GDCprepare(query_lgg_clin)

# Reshape proteins data (samples x (1+protein); first column id)
protein <- as_tibble(cbind(barcode = names(lgg_data)[-c(1:4)], 
                           t(lgg_data[,-c(1:4)])))
protein <- protein %>%
  janitor::row_to_names(1) %>% 
  mutate_at(2:nrow(protein), as.numeric) %>% 
  mutate_at(1, as.character) %>%
  mutate_at(1, str_replace, "-01A", "")

colnames(protein)[1] <- "barcode"

# Reshape clinical data 
clinical_drug <- lgg_data_clin$clinical_drug_lgg[-c(1:2),] %>%
  select(bcr_patient_barcode, bcr_drug_barcode, form_completion_date, 
         treatment_best_response, pharmaceutical_therapy_type) #%>%
  #filter(treatment_best_response != "[Not Available]") %>%
  #filter(treatment_best_response != "[Not Applicable]") %>%
  #filter(treatment_best_response != "[Unknown]")


clinical_drug <- clinical_drug %>%
  select(bcr_patient_barcode, treatment_best_response, pharmaceutical_therapy_type)#, treatment)

#clinical_drug <- clinical_drug[as.logical(1-duplicated(clinical_drug$bcr_patient_barcode)),]

clinical_patient <- lgg_data_clin$clinical_patient_lgg[-c(1:2),] %>%
  select(bcr_patient_barcode, treatment_outcome_first_course, gender, birth_days_to, race, #ethnicity,
         radiation_treatment_adjuvant)

clinical <- full_join(clinical_drug, clinical_patient, by = "bcr_patient_barcode")

#clinical <- clinical %>% 
#  filter(treatment_best_response != "[Not Available]") %>% 
#  filter(treatment_best_response != "[Not Applicable]") %>% 
#  filter(treatment_best_response != "[Unknown]")

colnames(clinical)[1] <- "barcode"

dat <- inner_join(clinical, protein, by = "barcode")
#dat <- full_join(clinical, protein, by = "barcode")

dat <- dat %>%
  filter(pharmaceutical_therapy_type != "[Not Available]") %>%
  filter(treatment_best_response != "[Not Available]") #%>%
  #filter(treatment_outcome_first_course != "[Not Available]") %>%

table(dat$treatment, dat$radiation_treatment_adjuvant)
table(dat$treatment)
for(i in 1:nrow(dat)){
  if((dat$radiation_treatment_adjuvant[i] == "YES") & (dat$treatment[i] == "Chemiotherapy")){
    dat$treatment[i] <- "Advanced"
  }
}
table(dat$treatment)

dat$treatment <- dat$treatment %>%
  replace(.=="Chemioterapy", "Standard") %>%
  replace(.=="Targeted Molecular therapy", "Advanced")

dat$treatment[67] <- "Advanced"#dat$radiation_treatment_adjuvant[67]

#samples_metadata <- TCGAbiolinks:::colDataPrepare(colnames(lgg_data[,-c(1:5)]))

#for(i in 1:435){
#  print(samples_metadata$treatments[[i]]$treatment_type)
#}

# se i trattamenti sono somminisrtati tutti lo stesso giorno e sono strategie 
# diverse (eg Targeted Molecular therapy and Chemotherapy) lo considero come 
# Targeted Molecular therapy

n <- length(unique(clinical_drug$bcr_patient_barcode))
therapy <- matrix(nrow=n, ncol=2)
for(i in 1:n){
  pat <- unique(clinical_drug$bcr_patient_barcode)[i]
  wt <- clinical_drug %>%
    filter(bcr_patient_barcode == pat)
  if(nrow(wt) > 1){
    if(length(unique(wt$pharmaceutical_therapy_type)) > 1){
      therapy[i,] <- c(pat, "Targeted Molecular therapy")
    } else {
      therapy[i,] <- c(pat, unique(wt$pharmaceutical_therapy_type))
    }
  } else {
    therapy[i,] <- c(pat, wt$pharmaceutical_therapy_type)
  }
}

therapy <- as_tibble(therapy)
colnames(therapy) <- c("bcr_patient_barcode", "treatment")
# ho rimosso i duplicati considerando solo quello "di sintesi"
clinical_drug <- clinical_drug[as.logical(1-duplicated(clinical_drug$bcr_patient_barcode)),]

clinical_drug <- full_join(clinical_drug, therapy, by = "bcr_patient_barcode")

