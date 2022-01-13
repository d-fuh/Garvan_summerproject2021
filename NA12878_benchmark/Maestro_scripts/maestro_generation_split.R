# This script generates a "maestro" file containing all annotated information from NA12878 SV callers and their subsets. 
# (Designated "small", "medium", and "large" by Coverage)

# NB: The final output has NOT been filtered by Quality. Please use the main workflow for manual filtering.
# Annotation_mode == "split". See more in ReadME.md.
# Final output: benchmark.split.maestro (data frame)
# START WORKFLOW

library(tidyverse)
library(vcfR)

#################################################################################################################################
################## Generate benchmark.split.maestro from AnnotSV-annotated files from all callers, split mode ###################
#################################################################################################################################


# reading in annotated TSV files
delly.split <- read_tsv("AnnotSV_NA12878_delly_split.tsv")
manta.split <- read_tsv("AnnotSV_NA12878_manta_split.tsv")
melt.split <- read_tsv("AnnotSV_NA12878_melt_split.tsv")
wham.split <- read_tsv("AnnotSV_NA12878_wham_split.tsv")
svaba.split <- read_tsv("AnnotSV_NA12878_svaba_split.tsv")

# Assigning callers
delly.split$Caller <- "Delly"
manta.split$Caller <- "Manta"
melt.split$Caller <- "Melt"
wham.split$Caller <- "Wham"
svaba.split$Caller <- "SvABA"


# joining all annotated files into one master file
library(plyr)

join(delly.split, manta.split, type = "full") %>%
  join(melt.split, type = "full") %>%
  join(wham.split, type = "full") %>%
  join(svaba.split, type = "full") -> benchmark.split 

# untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)


# optional clean-up to save memory
rm("delly.split"); rm("manta.split"); rm("melt.split"); rm("wham.split"); rm("svaba.split")


# append coverage data set information
benchmark.split$Coverage <- "Large"

# final filtering step
filter(benchmark.split, FILTER == "PASS") -> benchmark.split

# write resultant file into CSV
write.csv(benchmark.split, "maestro.split.csv")

# END WORKFLOW

