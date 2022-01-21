library(tidyverse)

###################################################################################################################
######################## Generate ExRes.csv from AnnotSV-annotated files from all samples #########################
###################################################################################################################
# Replace "ERXXX_XXXXX" with correct sample ID
# Sample list
## ER019_SAR1N
## ER052_MEL1N
## ER057_MEL2N
## ER095_MEL3N
## ER108_MEL5N
## ER118_MEL6N


# Reading in annotated files for each of the six ExRes samples
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")  ## sample 1
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")  ## sample 2, etc
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")
ERXXX_XXXXX <- read_tsv("ERXXX_XXXXX_annot.tsv")



# Joining all ExRes data into one master file
library(plyr)

join(ERXXX_XXXXX, ERXXX_XXXXX, type = "full") %>%  ## samples 1,2
  join(ERXXX_XXXXX, type = "full") %>%  ## samples 12, 3, etc
  join(ERXXX_XXXXX, type = "full") %>%  ## 123, 4
  join(ERXXX_XXXXX, type = "full") %>%  ## 1234, 5
  join(ERXXX_XXXXX, type = "full") -> ExRes.all

# Untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)

# Final filtering
filter(ExRes.all, FILTER == "PASS") -> ExRes.all  ## this step should be redundant, as each ExRes file came filtered.

# Write in csv for later use
write.csv(ExRes.all, "ExRes.csv")

# Clear intermediate files
rm(ExRes.all)
rm("ERXXX_XXXXX"); rm("ERXXX_XXXXX"); rm("ERXXX_XXXXX"); rm("ERXXX_XXXXX"); rm("ERXXX_XXXXX"); rm("ERXXX_XXXXX")
