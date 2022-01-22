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
ER019_SAR1N <- read_tsv("ER019_SAR1N_annot.tsv")  ## sample 1
ER052_MEL1N <- read_tsv("ER052_MEL1N_annot.tsv")  ## sample 2, etc
ER057_MEL2N <- read_tsv("ER057_MEL2N_annot.tsv")
ER095_MEL3N <- read_tsv("ER095_MEL3N_annot.tsv")
ER108_MEL5N <- read_tsv("ER108_MEL5N_annot.tsv")
ER118_MEL6N <- read_tsv("ER118_MEL6N_annot.tsv")


# Append sample ID to each frame before joining
ER019_SAR1N$ExResID <- "ER019_SAR1N"
ER052_MEL1N$ExResID <- "ER052_MEL1N"
ER057_MEL2N$ExResID <- "ER057_MEL2N"
ER095_MEL3N$ExResID <- "ER095_MEL3N"
ER108_MEL5N$ExResID <- "ER108_MEL5N"
ER118_MEL6N$ExResID <- "ER118_MEL6N"


# Joining all ExRes data into one master file
library(plyr)

join(ER019_SAR1N, ER052_MEL1N, type = "full") %>%  ## samples 1,2
  join(ER057_MEL2N, type = "full") %>%  ## samples 12, 3, etc
  join(ER095_MEL3N, type = "full") %>%  ## 123, 4
  join(ER108_MEL5N, type = "full") %>%  ## 1234, 5
  join(ER118_MEL6N, type = "full") -> ExRes.all

# Untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)

# Final filtering
filter(ExRes.all, FILTER == "PASS") -> ExRes.all  ## this step should be redundant, as each ExRes file came filtered.

# Write in csv for later use
write.csv(ExRes.all, "ExRes.csv")

# Clear intermediate files
rm(ExRes.all)
rm("ER019_SAR1N"); rm("ER052_MEL1N"); rm("ER057_MEL2N"); rm("ER095_MEL3N"); rm("ER108_MEL5N"); rm("ER118_MEL6N")
