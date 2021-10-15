# study-treatppmx
This repo contains all the scripts and the study needed for treatppmx project

## Structure of the project folder

### study-treatppmx
top-level folder and contains all of the folders and files associated with the project. 
### data
data contains the different scenarios used in the project. These files should not be altered and are ideally read-only. The actual processing of raw data is done in `treatppmx` package.

### doc
doc contains old scripts (simulation studies, theorethical studies etc) that I don't use anymore, but are still useful as guideline/reference.

### figs
figs contains any plots, images, tables, or figures created and saved by your code. It should be possible to delete and regenerate this folder with the scripts in the project folder.

### output
output contains several folder. Each folder stores the results (in `.rda` or `.RData` format) of a specific simulation scenario, sensitivity/parameter tuning study, theorethical studies. 

### src
src is a folder mainly containing files that are `source()`d in scripts stored in `doc`. 

#### more (for me)
https://kdestasio.github.io/post/r_best_practices/

cerca di creare file `info-nomecartella.md` per descrivere il contenuto di ogni folder/subfolder.
nelle cartelle dei risultati metti anche un file con la descrizione dell'ambiente R caricato