# This is the script version of the notebook
# "NA12878_benchmark.Rmd"
# for the 50% reads sample of NA12878

library(tidyverse)
library(GenomicRanges)
library(StructuralVariantAnnotation)

setwd("~/GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF")

# truth vcf
personalis_vcf <- readVcf("test_personalis_hg38.vcf")
spiral_vcf <- readVcf("test_spiral_hg38.vcf")

truth_vcf <- rbind(personalis_vcf, spiral_vcf)
truth_svgr <- breakpointRanges(truth_vcf)

# 50pc set
# original Dropbox file were compressed
# below files were decompressed
delly_vcf <- readVcf("NA12878_50pc_sampled.delly.vcf", "hg38")
manta_vcf <- readVcf("NA12878_50pc_sampled.manta.vcf", "hg38")
melt_vcf <- readVcf("NA12878_50pc_sampled.melt.vcf", "hg38")
wham_vcf <- readVcf("NA12878_50pc_sampled.wham.vcf", "hg38")
svaba_vcf <- readVcf("NA12878_50pc_sampled.svaba.vcf", "hg38") ## processed, converted file in folder


# bpr or svgr
delly_svgr <- breakpointRanges(delly_vcf)
manta_svgr <- breakpointRanges(manta_vcf)
melt_svgr <- breakpointRanges(melt_vcf)
wham_svgr <- breakpointRanges(wham_vcf)
svaba_svgr <- breakpointRanges(svaba_vcf)

# filter
delly_svgr=subset(delly_svgr, delly_svgr@elementMetadata[, 5] == "PASS") 
manta_svgr=subset(manta_svgr, manta_svgr@elementMetadata[, 5] == "PASS") 
melt_svgr=subset(melt_svgr, melt_svgr@elementMetadata[, 5] == "PASS") 
wham_svgr=subset(wham_svgr, wham_svgr@elementMetadata[, 5] == "PASS") 
svaba_svgr=subset(svaba_svgr, svaba_svgr@elementMetadata[, 5] == "PASS") 

# assign caller info
delly_svgr$Caller <- "Delly"
manta_svgr$Caller <- "Manta"
melt_svgr$Caller <- "Melt"
wham_svgr$Caller <- "Wham"
svaba_svgr$Caller <- "SvABA"

# assign coverage info
delly_svgr$Cov <- "50pc"
manta_svgr$Cov <- "50pc"
melt_svgr$Cov <- "50pc"
wham_svgr$Cov <- "50pc"
svaba_svgr$Cov <- "50pc"

# main svgr
svgr <- c(delly_svgr, manta_svgr, melt_svgr, wham_svgr, svaba_svgr)

svgr$truth_matches <- countBreakpointOverlaps(svgr, truth_svgr, 
                                              
                                              maxgap=100, sizemargin=0.25, ## explain
                                              
                                              restrictMarginToSizeMultiple=0.5, ## explain
                                              
                                              countOnlyBest=TRUE)
# get summary stats
fiftypc = as.data.frame(svgr) %>%
  dplyr::select(Caller, truth_matches) %>%
  dplyr::group_by(Caller) %>%
  dplyr::summarise(
    calls=dplyr::n(), ## number of calls for each caller that matches the truth
    tp=sum(truth_matches > 0)) %>% ## tp=true positive calls for each caller = sum(all non-0 tp calls for each caller)
  
  dplyr::group_by(Caller) %>%
  
  dplyr::mutate(
    cum_tp=cumsum(tp),
    cum_n=cumsum(calls),
    cum_fp=cum_n - cum_tp,
    
    Precision=cum_tp / cum_n,
    Recall=cum_tp/length(truth_svgr)
  )
fiftypc$Cov <- "50pc"

# plot
ggplot(fiftypc) +
  geom_point(aes(x=Caller, y=Recall), shape=25, fill="red", size=2.5) +
  geom_point(aes(x=Caller, y=Precision), shape=3, size=3) +
  facet_wrap(~Cov) +
  ylab("Precision (cross) & Recall (Red triangle)")








