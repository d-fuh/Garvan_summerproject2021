# This script generates a "maestro" file containing annotated information from NA12878 SV callers and their subsets
# using the online tool AnnotSV. 
# (Files and respective data frames are designated "small", "medium", and "large" by Coverage)

# Annotation_mode == "full". See more in ReadME.md.
# Final output: {benchmark.maestro (data frame), maestro.full.csv}
# START WORKFLOW

library(tidyverse)

#################################################################################################################################
######################## Generate benchmark.annot.maestro from AnnotSV-annotated files from all callers #########################
#################################################################################################################################

# reading in annotated TSV files
## small
s.delly.annot <- read_tsv("AnnotSV_NA12878_small_delly.tsv")
s.manta.annot <- read_tsv("AnnotSV_NA12878_small_manta.tsv")
s.melt.annot <- read_tsv("AnnotSV_NA12878_small_melt.tsv")
s.wham.annot <- read_tsv("AnnotSV_NA12878_small_wham.tsv")
s.svaba.annot <- read_tsv("AnnotSV_NA12878_small_svaba.tsv")
## medium
m.delly.annot <- read_tsv("AnnotSV_NA12878_med_delly.tsv")
m.manta.annot <- read_tsv("AnnotSV_NA12878_med_manta.tsv")
m.melt.annot <- read_tsv("AnnotSV_NA12878_med_melt.tsv")
m.wham.annot <- read_tsv("AnnotSV_NA12878_med_wham.tsv")
m.svaba.annot <- read_tsv("AnnotSV_NA12878_med_svaba.tsv")
## large
delly.annot <- read_tsv("AnnotSV_NA12878_delly.tsv")
manta.annot <- read_tsv("AnnotSV_NA12878_manta.tsv")
melt.annot <- read_tsv("AnnotSV_NA12878_melt.tsv")
wham.annot <- read_tsv("AnnotSV_NA12878_wham.tsv")
svaba.annot <- read_tsv("AnnotSV_NA12878_svaba.tsv")

# Assigning callers
s.delly.annot$Caller <- "Delly"
s.manta.annot$Caller <- "Manta"
s.melt.annot$Caller <- "Melt"
s.wham.annot$Caller <- "Wham"
s.svaba.annot$Caller <- "SvABA"

m.delly.annot$Caller <- "Delly"
m.manta.annot$Caller <- "Manta"
m.melt.annot$Caller <- "Melt"
m.wham.annot$Caller <- "Wham"
m.svaba.annot$Caller <- "SvABA"

delly.annot$Caller <- "Delly"
manta.annot$Caller <- "Manta"
melt.annot$Caller <- "Melt"
wham.annot$Caller <- "Wham"
svaba.annot$Caller <- "SvABA"


# joining all annotated files into one master file
library(plyr)
## small
join(s.delly.annot, s.manta.annot, type = "full") %>%
  join(s.melt.annot, type = "full") %>%
  join(s.wham.annot, type = "full") %>%
  join(s.svaba.annot, type = "full") -> s.benchmark.annot

## medium
join(m.delly.annot, m.manta.annot, type = "full") %>%
  join(m.melt.annot, type = "full") %>%
  join(m.wham.annot, type = "full") %>%
  join(m.svaba.annot, type = "full") -> m.benchmark.annot

## large: r.d.s.
join(delly.annot, manta.annot, type = "full") %>%
  join(melt.annot, type = "full") %>%
  join(wham.annot, type = "full") %>%
  join(svaba.annot, type = "full") -> benchmark.annot 

# untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)


# optional clean-up to save memory

rm("s.delly.annot"); rm("s.manta.annot"); rm("s.melt.annot"); rm("s.wham.annot"); rm("s.svaba.annot")

rm("m.delly.annot"); rm("m.manta.annot"); rm("m.melt.annot"); rm("m.wham.annot"); rm("m.svaba.annot")

rm("delly.annot"); rm("manta.annot"); rm("melt.annot"); rm("wham.annot"); rm("svaba.annot")


# append coverage dataset information
s.benchmark.annot$Coverage <- "Small"
m.benchmark.annot$Coverage <- "Medium"
benchmark.annot$Coverage <- "Large"

# final filtering step
filter(s.benchmark.annot, FILTER == "PASS") -> s.benchmark.annot
filter(m.benchmark.annot, FILTER == "PASS") -> m.benchmark.annot
filter(benchmark.annot, FILTER == "PASS") -> benchmark.annot


# combine Annot files from each coverage set (small, medium, large) into one master file

library(plyr)

join(s.benchmark.annot, m.benchmark.annot, type = "full") %>%
  join(benchmark.annot, type = "full") -> benchmark.full

detach("package:plyr", unload = TRUE)

write.csv(benchmark.full, "maestro.full.csv")

# clear intermediate files
rm("s.benchmark.annot"); rm("m.benchmark.annot"); rm("benchmark.annot")
