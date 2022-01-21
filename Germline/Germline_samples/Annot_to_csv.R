library(tidyverse)

###################################################################################################################
######################## Generate ERXXX.csv from AnnotSV-annotated files from all callers #########################
###################################################################################################################
# Replace "ERXXX_XXXXX" with correct sample ID

# Reading in annotated files
delly.annot <- read_tsv("AnnotSV_ERXXX_XXXXX_delly.tsv")
manta.annot <- read_tsv("AnnotSV_ERXXX_XXXXX_manta.tsv")
melt.annot <- read_tsv("AnnotSV_ERXXX_XXXXX_melt.tsv")
wham.annot <- read_tsv("AnnotSV_ERXXX_XXXXX_wham.tsv")
svaba.annot <- read_tsv("AnnotSV_ERXXX_XXXXX_svaba.tsv")

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
filter(all.annot, FILTER == "PASS") -> ERXXX_XXXXX_annot

# Write in csv for later use
write.csv(ERXXX_XXXXX_annot, "ERXXX_XXXXX_annot.csv")

# Clear intermediate files
rm(all.annot)
rm(ERXXX_XXXXX_annot)
