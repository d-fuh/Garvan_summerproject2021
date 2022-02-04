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
ft_delly_vcf <- readVcf("NA12878_50pc_sampled.delly.vcf", "hg38")
ft_manta_vcf <- readVcf("NA12878_50pc_sampled.manta.vcf", "hg38")
ft_melt_vcf <- readVcf("NA12878_50pc_sampled.melt.vcf", "hg38")
ft_wham_vcf <- readVcf("NA12878_50pc_sampled.wham.vcf", "hg38")
ft_svaba_vcf <- readVcf("NA12878_50pc_sampled.svaba.vcf", "hg38") ## processed, converted file in folder


# bpr or svgr
ft_delly_svgr <- breakpointRanges(ft_delly_vcf)
ft_manta_svgr <- breakpointRanges(ft_manta_vcf)
ft_melt_svgr <- breakpointRanges(ft_melt_vcf)
ft_wham_svgr <- breakpointRanges(ft_wham_vcf)
ft_svaba_svgr <- breakpointRanges(ft_svaba_vcf)

# filter
ft_delly_svgr=subset(ft_delly_svgr, ft_delly_svgr@elementMetadata[, 5] == "PASS") 
ft_manta_svgr=subset(ft_manta_svgr, ft_manta_svgr@elementMetadata[, 5] == "PASS") 
ft_melt_svgr=subset(ft_melt_svgr, ft_melt_svgr@elementMetadata[, 5] == "PASS") 
ft_wham_svgr=subset(ft_wham_svgr, ft_wham_svgr@elementMetadata[, 5] == "PASS") 
ft_svaba_svgr=subset(ft_svaba_svgr, ft_svaba_svgr@elementMetadata[, 5] == "PASS") 

# assign caller info
ft_delly_svgr$Caller <- "Delly"
ft_manta_svgr$Caller <- "Manta"
ft_melt_svgr$Caller <- "Melt"
ft_wham_svgr$Caller <- "Wham"
ft_svaba_svgr$Caller <- "SvABA"

# assign coverage info
ft_delly_svgr$Cov <- "18x"
ft_manta_svgr$Cov <- "18x"
ft_melt_svgr$Cov <- "18x"
ft_wham_svgr$Cov <- "18x"
ft_svaba_svgr$Cov <- "18x"




# truth matches
ft_delly_svgr$truth_matches <- countBreakpointOverlaps(ft_delly_svgr, truth_svgr, maxgap=200, sizemargin=0.25,
                                                          restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
ft_manta_svgr$truth_matches <- countBreakpointOverlaps(ft_manta_svgr, truth_svgr, maxgap=200, sizemargin=0.25,
                                                          restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
ft_melt_svgr$truth_matches <- countBreakpointOverlaps(ft_melt_svgr, truth_svgr, maxgap=200, sizemargin=0.25,
                                                         restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
ft_wham_svgr$truth_matches <- countBreakpointOverlaps(ft_wham_svgr, truth_svgr, maxgap=200, sizemargin=0.25,
                                                         restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
ft_svaba_svgr$truth_matches <- countBreakpointOverlaps(ft_svaba_svgr, truth_svgr, maxgap=200, sizemargin=0.25,
                                                          restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
# main svgr
ft_svgr <- c(ft_delly_svgr, ft_manta_svgr, ft_melt_svgr, ft_wham_svgr, ft_svaba_svgr)

# get summary stats
PR_fiftypc = as.data.frame(ft_svgr) %>%
  dplyr::select(Caller, truth_matches) %>%
  dplyr::group_by(Caller) %>%
  dplyr::summarise(
    calls=dplyr::n(), ## number of calls for each caller that matches the truth
    tp=sum(truth_matches > 0)) %>% ## tp=true positive calls for each caller = sum(all non-0 tp calls for each caller)
  
  dplyr::mutate(
    fp=calls-tp,
    Precision=tp/calls,
    Recall=tp/length(truth_svgr)
  )

PR_fiftypc$Cov <- "18x"

# plot
ggplot(PR_fiftypc) +
  geom_point(aes(x=Caller, y=Recall), shape=25, fill="red", size=2.5) +
  geom_point(aes(x=Caller, y=Precision), shape=3, size=3) +
  facet_wrap(~Cov) +
  ylab("Precision (cross) & Recall (Red triangle)")





