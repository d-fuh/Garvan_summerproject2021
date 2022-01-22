library(tidyverse)

###################################################################################################################
######################## Generate ERXXX.csv from AnnotSV-annotated files from all callers #########################
###################################################################################################################
# Replace "ER118_MEL6N" with correct sample ID
# Sample list
## ER019_SAR1N [v]
## ER052_MEL1N [v]
## ER057_MEL2N [v]
## ER095_MEL3N [v]
## ER108_MEL5N [v]
## ER118_MEL6N [v]

# Reading in annotated files
delly.annot <- read_tsv("AnnotSV_ER118_MEL6N_delly.tsv")
manta.annot <- read_tsv("AnnotSV_ER118_MEL6N_manta.tsv")
melt.annot <- read_tsv("AnnotSV_ER118_MEL6N_melt.tsv")
wham.annot <- read_tsv("AnnotSV_ER118_MEL6N_wham.tsv")
svaba.annot <- read_tsv("AnnotSV_ER118_MEL6N_svaba.tsv")

# Assigning callers
delly.annot$Caller <- "Delly"
manta.annot$Caller <- "Manta"
melt.annot$Caller <- "Melt"
wham.annot$Caller <- "Wham"
svaba.annot$Caller <- "SvABA"


# Joining all annotated files into one master file
library(plyr)

join(delly.annot, manta.annot, type = "full") %>%
  join(melt.annot, type = "full") %>%
  join(wham.annot, type = "full") %>%
  join(svaba.annot, type = "full") -> all.annot

# Untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)


# Clean-up to save memory
rm("delly.annot"); rm("manta.annot"); rm("melt.annot"); rm("wham.annot"); rm("svaba.annot")


# Final filtering
filter(all.annot, FILTER == "PASS") -> ER118_MEL6N_annot

# Write in csv for later use
write.csv(ER118_MEL6N_annot, "ER118_MEL6N_annot.csv")

# Clear intermediate files
rm(all.annot)
rm(ER118_MEL6N_annot)
