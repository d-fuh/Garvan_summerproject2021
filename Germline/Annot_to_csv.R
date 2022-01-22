library(tidyverse)

###################################################################################################################
######################## Generate ERXXX.csv from AnnotSV-annotated files from all callers #########################
###################################################################################################################
# Replace "ER019_SAR1N" with correct sample ID
# Sample list
## ER019_SAR1N
## ER052_MEL1N
## ER057_MEL2N
## ER095_MEL3N
## ER108_MEL5N
## ER118_MEL6N

# Reading in annotated files
delly.annot <- read_tsv("AnnotSV_ER019_SAR1N_delly.tsv")
manta.annot <- read_tsv("AnnotSV_ER019_SAR1N_manta.tsv")
melt.annot <- read_tsv("AnnotSV_ER019_SAR1N_melt.tsv")
wham.annot <- read_tsv("AnnotSV_ER019_SAR1N_wham.tsv")
svaba.annot <- read_tsv("AnnotSV_ER019_SAR1N_svaba.tsv")

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
filter(all.annot, FILTER == "PASS") -> ER019_SAR1N_annot

# Write in csv for later use
write.csv(ER019_SAR1N_annot, "ER019_SAR1N_annot.csv")

# Clear intermediate files
rm(all.annot)
rm(ER019_SAR1N_annot)
