# In this script, we will output a benchmark result from NA12878 small and medium subset using our previous workflow.
# after which, we will combine these results in order to create a wrapped bar plot for counts ~ SV_type, group_by = coverage.

# NB: The final output "benchmark.annot.maestro" has NOT been filtered by quality. Please use the main workflow for manual filtering.
# START WORKFLOW

library(tidyverse)
library(vcfR)

# load VCF
## small
s.delly.vcf <- read.vcfR("NA12878_small.delly.vcf")
s.manta.vcf <- read.vcfR("NA12878_small.manta.vcf")
s.melt.vcf <- read.vcfR("NA12878_small.melt.vcf")
s.wham.vcf <- read.vcfR("NA12878_small.wham.vcf") 
s.svaba.vcf <- read.vcfR("NA12878_small.svaba.vcf")

## medium
m.delly.vcf <- read.vcfR("NA12878_med.delly.vcf")
m.manta.vcf <- read.vcfR("NA12878_med.manta.vcf")
m.melt.vcf <- read.vcfR("NA12878_med.melt.vcf")
m.wham.vcf <- read.vcfR("NA12878_med.wham.vcf") 
m.svaba.vcf <- read.vcfR("NA12878_med.svaba.vcf")

## large
delly.vcf <- read.vcfR("NA12878.delly.vcf")
manta.vcf <- read.vcfR("NA12878.manta.vcf")
melt.vcf <- read.vcfR("NA12878.melt.vcf")
wham.vcf <- read.vcfR("NA12878.wham.vcf") 
svaba.vcf <- read.vcfR("NA12878.svaba.vcf")


# filtering only PASS SVs
## NB: Is "getFIX" correct?
## small
as.data.frame(getFIX(s.delly.vcf)) %>% 
  filter(FILTER == "PASS") -> s.delly.pass
as.data.frame(getFIX(s.manta.vcf)) %>% 
  filter(FILTER == "PASS") -> s.manta.pass
as.data.frame(getFIX(s.melt.vcf)) %>% 
  filter(FILTER == "PASS") -> s.melt.pass
as.data.frame(getFIX(s.wham.vcf)) %>% 
  filter(FILTER == "PASS") -> s.wham.pass
as.data.frame(getFIX(s.svaba.vcf)) %>% 
  filter(FILTER == "PASS") -> s.svaba.pass

## medium
as.data.frame(getFIX(m.delly.vcf)) %>% 
  filter(FILTER == "PASS") -> m.delly.pass
as.data.frame(getFIX(m.manta.vcf)) %>% 
  filter(FILTER == "PASS") -> m.manta.pass
as.data.frame(getFIX(m.melt.vcf)) %>% 
  filter(FILTER == "PASS") -> m.melt.pass
as.data.frame(getFIX(m.wham.vcf)) %>% 
  filter(FILTER == "PASS") -> m.wham.pass
as.data.frame(getFIX(m.svaba.vcf)) %>% 
  filter(FILTER == "PASS") -> m.svaba.pass

## large
as.data.frame(getFIX(delly.vcf)) %>% 
  filter(FILTER == "PASS") -> delly.pass
as.data.frame(getFIX(manta.vcf)) %>% 
  filter(FILTER == "PASS") -> manta.pass
as.data.frame(getFIX(melt.vcf)) %>% 
  filter(FILTER == "PASS") -> melt.pass
as.data.frame(getFIX(wham.vcf)) %>% 
  filter(FILTER == "PASS") -> wham.pass
as.data.frame(getFIX(svaba.vcf)) %>% 
  filter(FILTER == "PASS") -> svaba.pass

# clear
rm("s.delly.vcf"); rm("s.manta.vcf"); rm("s.melt.vcf"); rm("s.wham.vcf"); rm("s.svaba.vcf")

rm("m.delly.vcf"); rm("m.manta.vcf"); rm("m.melt.vcf"); rm("m.wham.vcf"); rm("m.svaba.vcf")

rm("delly.vcf"); rm("manta.vcf"); rm("melt.vcf"); rm("wham.vcf"); rm("svaba.vcf")


# assigning respective callers
## small
s.delly.pass$Caller <- "Delly"
s.manta.pass$Caller <- "Manta"
s.melt.pass$Caller <- "Melt"
s.wham.pass$Caller <- "Wham"
s.svaba.pass$Caller <- "SvABA"

## medium
m.delly.pass$Caller <- "Delly"
m.manta.pass$Caller <- "Manta"
m.melt.pass$Caller <- "Melt"
m.wham.pass$Caller <- "Wham"
m.svaba.pass$Caller <- "SvABA"

## large
delly.pass$Caller <- "Delly"
manta.pass$Caller <- "Manta"
melt.pass$Caller <- "Melt"
wham.pass$Caller <- "Wham"
svaba.pass$Caller <- "SvABA"

# synthesise the filtered data frames into one master file
library(plyr)
## small
join(s.delly.pass, s.manta.pass, type = "full") %>%
  join(s.melt.pass, type = "full") %>%
  join(s.wham.pass, type = "full") %>%
  join(s.svaba.pass, type = "full") -> s.benchmark.master

## medium
join(m.delly.pass, m.manta.pass, type = "full") %>%
  join(m.melt.pass, type = "full") %>%
  join(m.wham.pass, type = "full") %>%
  join(m.svaba.pass, type = "full") -> m.benchmark.master

## large
join(delly.pass, manta.pass, type = "full") %>%
  join(melt.pass, type = "full") %>%
  join(wham.pass, type = "full") %>%
  join(svaba.pass, type = "full") -> benchmark.master

# untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)

rm("s.delly.pass"); rm("s.manta.pass"); rm("s.melt.pass"); rm("s.wham.pass"); rm("s.svaba.pass")

rm("m.delly.pass"); rm("m.manta.pass"); rm("m.melt.pass"); rm("m.wham.pass"); rm("m.svaba.pass")

rm("delly.pass"); rm("manta.pass"); rm("melt.pass"); rm("wham.pass"); rm("svaba.pass")


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

## large
join(delly.annot, manta.annot, type = "full") %>%
  join(melt.annot, type = "full") %>%
  join(wham.annot, type = "full") %>%
  join(svaba.annot, type = "full") -> benchmark.annot 

##TODO double check that this is the correct way to join these data frames....

# untick plyr afterwards to avoid conflict with dplyr
detach("package:plyr", unload = TRUE)


# optional clean-up to save memory

rm("s.delly.annot"); rm("s.manta.annot"); rm("s.melt.annot"); rm("s.wham.annot"); rm("s.svaba.annot")

rm("m.delly.annot"); rm("m.manta.annot"); rm("m.melt.annot"); rm("m.wham.annot"); rm("m.svaba.annot")

rm("delly.annot"); rm("manta.annot"); rm("melt.annot"); rm("wham.annot"); rm("svaba.annot")

# combine annot files from each coverage set (small, medium, large) into one master file

library(plyr)

join(s.benchmark.annot, m.benchmark.annot, type = "full") %>%
  join(benchmark.annot, type = "full") -> benchmark.annot.maestro

detach("package:plyr", unload = TRUE)

write.csv(benchmark.annot.maestro, "maestro.annot.csv")

# OPTIONAL: USE IF NOT USING NA12878_benchmark.rmd 
# rm("s.benchmark.annot"); rm("m.benchmark.annot"); rm("benchmark.annot")
# rm("s.benchmark.master"); rm("m.benchmark.master"); rm("benchmark.master")
