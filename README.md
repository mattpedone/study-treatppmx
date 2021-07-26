# study-treatppmx
This repo contains all the scripts and the study needed for treatppmx project

## Structure of the project folder

### study-treatppmx
top-level folder and contains all of the folders and files associated with the project. 

### data
data contains the raw data files used in the project. These files should not be altered and are ideally read-only. Some may be duplicated as present also in treatppmx package

### doc
doc contains any manuscripts or interim summaries produced with the project. it is going to contain several folders (simulation studies, theorethical studies etc)

### figs
figs contains any plots, images, tables, or figures created and saved by your code. It should be possible to delete and regenerate this folder with the scripts in the project folder.

### output
output contains non-figure objects created by the scripts. For example, processed data or logs.

### src
src is an optional folder for any files you may want to source() in your scripts. This is not code that is run. For example, simple .R files containing functions.

#### more (for me)
https://kdestasio.github.io/post/r_best_practices/
