# Generate GRange object from VCFs

library(tidyverse)
library(vcfR)
library(GenomicRanges)

##################################################################################################################################
######################################################### DATA PREP ##############################################################
##################################################################################################################################

# Can we generate gr from annot.tsv?
## read in annotated SV call set from all callers from all subsets of NA12878 (WARNING: overlapping subsets)
benchmark.full <- read.csv("./Maestro_scripts/maestro.full.csv")

## filtering overlapping-subset master file to purer call set with all adequate info that we will use from here on
large.full <- filter(benchmark.full, Coverage == "Large") # careful: cap sensitive

## further subsetting the call set to extract infos we need
large.tst <- large.full %>% 
  subset(., select = c("SV_chrom", "SV_start", "SV_end", "QUAL", "Caller", "SV_type"))

## always check dim
dim(benchmark.full); dim(large.full); dim(large.tst)

# generate gr from delly.vcf as an example
## delly
delly.tst <- filter(large.tst, Caller == "Delly")

delly.gr <- GRanges(
  seqnames = Rle(delly.tst$SV_chrom, delly.tst$rownames),
  ranges = IRanges(delly.tst$SV_start, end = delly.tst$SV_end, names = delly.tst$rownames),
)

## manta
manta.tst <- filter(large.tst, Caller == "Manta")

manta.gr <- GRanges(
  seqnames = Rle(manta.tst$SV_chrom, manta.tst$rownames),
  ranges = IRanges(manta.tst$SV_start, end = manta.tst$SV_end, names = manta.tst$rownames),
)

## svaba
svaba.tst <- filter(large.tst, Caller == "SvABA")

svaba.gr <- GRanges(
  seqnames = Rle(svaba.tst$SV_chrom, svaba.tst$rownames),
  ranges = IRanges(svaba.tst$SV_start, end = svaba.tst$SV_end, names = svaba.tst$rownames),
)

## melt
melt.tst <- filter(large.tst, Caller == "Melt")

melt.gr <- GRanges(
  seqnames = Rle(melt.tst$SV_chrom, melt.tst$rownames),
  ranges = IRanges(melt.tst$SV_start, end = melt.tst$SV_end, names = melt.tst$rownames),
)

## wham
wham.tst <- filter(large.tst, Caller == "Wham")

wham.gr <- GRanges(
  seqnames = Rle(wham.tst$SV_chrom, wham.tst$rownames),
  ranges = IRanges(wham.tst$SV_start, end = wham.tst$SV_end, names = wham.tst$rownames),
)


##################################################################################################################################
############################################################  MAIN  ##############################################################
##################################################################################################################################

# Define two matching SV calls as having "similar SV range". By "similar", we define their boundaries (SV_start & SV_end) to be
# within each other's mean SV length by 0.2%.

# Graphical example. e.g.

# 5'---------------*|------------------------A-----------------|*----------------------------------------3'
# 5'---------------------|--------------------B--------------------|*------------------------------------3'
# 5'---------------*|---------------------------------C----------------------------------|*--------------3'

# * = breakpoint. A, B, C are SV calls each defined by their start & end positions.

# Q: Can A, B, C be considered the same calls?

# A: We may require for two calls to be considered the same, they must have their SV_start = (another's SV_start) +- some limit
# , and similarly for their SV_end positions.
# The two conditions (that the two's start AND end are within certain bounds) must be met before the two calls are considered
# "the same".

# To formalise this setting: Out goal is then

# > For every two matching SV calls i, j, each defined by their (SV_start[x], SV_end[x]), they have:

# (1) SV_start[i] = SV_start[j] "loosely"; &&
# (2) SV_end[i] = SV_end[j] "loosely".

# i.e. 
# (1) between(SV_start[i], SV_start[j]-CI, SV_start[j]+CI) == TRUE ## prove that the converse holds true automatically
# &&
# (2) between(SV_end[i], SV_end[j]-CI, SV_end[j]+CI) == TRUE

# set CI = 0.002*mean(SV_length[i,j]) for example


# > For such two matching calls, we group them by their Caller profiles:

# if((1) == TRUE && (2) == TRUE) {
#   aggregate(SV[i,j], by = Caller, collapse = ',')
# }



# find overlaps
findOverlaps(delly.gr, manta.gr) ## order: query, subject

