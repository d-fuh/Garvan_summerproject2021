# This is the script version of the notebook
# "NA12878_benchmark.Rmd"

library(tidyverse)
#library(paletteer)
#library(patchwork)
library(GenomicRanges)
library(StructuralVariantAnnotation)

setwd("~/GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF")

# truth vcf
personalis_vcf <- readVcf("test_personalis_hg38.vcf")
spiral_vcf <- readVcf("test_spiral_hg38.vcf")

truth_vcf <- rbind(personalis_vcf, spiral_vcf)

# 50pc set
delly_vcf <- readVcf("NA12878_50pc_sampled.delly.vcf", "hg38")
manta_vcf <- readVcf("NA12878_50pc_sampled.manta.vcf", "hg38")
melt_vcf <- readVcf("NA12878_50pc_sampled.melt.vcf", "hg38")
wham_vcf <- readVcf("NA12878_50pc_sampled.wham.vcf", "hg38")
svaba_vcf <- readVcf("NA12878_50pc_sampled.svaba.sv.vcf", "hg38") ## is this processed?


# bpr or svgr
delly_svgr <- breakpointRanges(delly_vcf)
manta_svgr <- breakpointRanges(manta_vcf)
melt_svgr <- breakpointRanges(melt_vcf)
wham_svgr <- breakpointRanges(wham_vcf)
svaba_svgr <- breakpointRanges(svaba_vcf)


