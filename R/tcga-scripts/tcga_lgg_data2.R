#BiocManager::install('GenomicDataCommons')
library(GenomicDataCommons)
GenomicDataCommons::status()
available_fields('projects') %>% length()

pQuery = projects()
default_fields(pQuery)

presults = pQuery %>% results_all() # or optionally results_all()
ids(presults) # numero 59

#Use available fields and grep fields to find the facets you want.
available_fields(files())
grep_fields('files','analysis.workflow')

#Use aggregrate to summarize the total number of files for each facet.
res = files() %>%
  #Default is to set facets from default_fields()
  facet(c('type','data_type','data_format',
          'cases.project.project_id')) %>%
  aggregations()

str(res)
res$cases.project.project_id$key  %>% .[order(.)]
res$data_type$key %>% .[order(.)]

res.projects <- projects() %>%
  facet(c("project_id")) %>%
  aggregations()

# str(res.projects)
res.projects

res.cases <- cases() %>%
  facet() %>% #Default is to set facets for all default fields.
  aggregations()

head(res.cases$primary_site)

desired_fields <-c("cases.project.project_id", default_fields('files'), 
                   grep_fields('files', "associated_entities"),  
                   "analysis.analysis_type", "analysis.workflow_type", 
                   "analysis.workflow_version")
#desired_fields <- c("cases.project.name", "cases.data_category")
length(desired_fields) #58 fields


qfiles <- files(fields=desired_fields) %>%
  filter(~ type == 'protein_expression' &
           #analysis.workflow_type == 'HTSeq - Counts' &
           (cases.project.project_id == "TCGA-LGG"))


qfiles %>% count()
str(qfiles)

res.expn <- results_all(x=qfiles)
head(res.expn$file_name)
head(res.expn$file_id)

# this checks if each file is associated with a single sample id
idx <- sapply(res.expn$associated_entities , nrow) %>% grep(2, .)

manifest_df = qfiles %>% manifest()
head(manifest_df)
dim(manifest_df)

write.table(manifest_df, "TARGET_Manifest_RNAseq_Counts.txt", row.names = FALSE, sep="\t", quote=FALSE)

options(gdc_client="/home/matt/Dropbox/PHD/study-treatppmx/gdc-client")
gdc_client()  
dir.create("Expn_Data")
gdc_set_cache(directory = "Expn_Data/")
fnames = gdcdata(manifest_df$id,progress=FALSE,access_method = "api", use_cached = FALSE)
head(fnames)

source("Cat_Expn_Data.r")
