# support read info extraction script
# for NA12878 small subset
library(StructuralVariantAnnotation)
setwd("~/GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF")


# Below are the stats we want to extract

# |Caller| SV support        | Total depth             | Note                 |
# |:-----|:------------------|:-----------------------:|:---------------------|
# |Delly | DR                | DR + DV                 |                      |
# |SvABA | AD                | DP                      | DP=total depth       |
# |Manta | PR[2] + SR[2]     | SR[1]+SR[2]+PR[1]+PR[2] | SR lists both the alt & ref so you need to split it  |
# |Wham  | A                 | n/a                     |                      |
# |Melt  | LP, RP            | LP + RP                 |                      |



## Delly: DR & (DR+DV)
delly_vcf<-readVcf("NA12878_small.delly.vcf")

# DR & DV are both in geno so import to INFO
#info(delly_vcf)$DR=geno(delly_vcf)$DR
#info(delly_vcf)$DV=geno(delly_vcf)$DV
  ##can check
  ##rownames(geno(delly_vcf)$DR)==rownames(info(delly_vcf))

# consistent naming is better
info(delly_vcf)$support=geno(delly_vcf)$DR
info(delly_vcf)$tot_depth=geno(delly_vcf)$DR + geno(delly_vcf)$DV

# create svgr object
delly.bpr=breakpointRanges(delly_vcf, info_columns=c("support", "tot_depth"))

# Assign a column for caller name
delly.bpr$Caller="Delly"

# Need to filter out the variants that pass
delly.bpr=delly.bpr[which(delly.bpr$FILTER=="PASS")]




## Svaba: AD & DP
svaba_vcf<-readVcf("NA12878_small.svaba.vcf")

# original svaba vcf has $SVTYPE == "BND" for all calls
# used a script to convert them into normal SV_type
# and written into svaba_converted.csv files
small.conv.svaba <- read.csv("small_svaba_converted.csv") %>% subset(select = -X)
info(svaba_vcf)$SVTYPE=small.conv.svaba$SV_type ## double check the ALT column matches
    #svaba_vcf@fixed$ALT == small.conv.svaba$ALT

# import support read info
#info(svaba_vcf)$AD=geno(svaba_vcf)$AD
#info(svaba_vcf)$DP=geno(svaba_vcf)$DP

# check rowname ID matches
  #rownames(AD) == rownames(info(svaba_vcf))

# consistent naming
info(svaba_vcf)$support=geno(svaba_vcf)$AD
info(svaba_vcf)$tot_depth=geno(svaba_vcf)$AD + geno(svaba_vcf)$DP

# create svgr object
svaba.bpr<-breakpointRanges(svaba_vcf, info_columns=c("support", "tot_depth"))

# Assign a column for caller name
svaba.bpr$Caller="SvABA"

# Need to filter out the variants that pass
svaba.bpr=svaba.bpr[which(svaba.bpr$FILTER=="PASS")]




## Manta: See summary at top
manta_vcf<-readVcf("NA12878_small.manta.vcf")
    #head(manta_vcf)

# PR, SR are in GT
# First, extract them
# PR
PR=geno(manta_vcf)$PR

## Note this entire thing is a list, You will need to extract every first item in this list.
info(manta_vcf)$PR_ref=sapply(PR, function(x) x[1]) ## PR[1]
info(manta_vcf)$PR_alt=sapply(PR, function(x) x[2]) ## PR[2]

# Same for SR
SR=geno(manta_vcf)$SR

## Note this entire thing is a list, You will need to extract every first item in this list.
info(manta_vcf)$SR_ref=sapply(SR, function(x) x[1]) ## SR[1]
info(manta_vcf)$SR_alt=sapply(SR, function(x) x[2]) ## SR[2]


# Then compute our stats
#
# SV support:   PR[2] + SR[2] 
# Total depth:  SR[1]+SR[2]+PR[1]+PR[2]


info(manta_vcf)$support=info(manta_vcf)$PR_alt+info(manta_vcf)$SR_alt
info(manta_vcf)$tot_depth=info(manta_vcf)$SR_ref+info(manta_vcf)$SR_alt+info(manta_vcf)$PR_ref+info(manta_vcf)$PR_alt


# create svgr object
  # omit BND_DEPTH
manta.bpr=breakpointRanges(manta_vcf, info_columns=c("support", "tot_depth"))


# assign caller
manta.bpr$Caller="Manta"

# filter
manta.bpr=manta.bpr[which(manta.bpr$FILTER=="PASS")]

rm(PR, SR)




## Wham: A
wham_vcf<-readVcf("NA12878_small.wham.vcf")

info(wham_vcf)$support=info(wham_vcf)$A

wham.bpr=breakpointRanges(wham_vcf, info_columns=c("support"))

# assign caller
wham.bpr$Caller="Wham"

# filter
wham.bpr=wham.bpr[which(wham.bpr$FILTER=="PASS")]




## Melt: LP, RP; LP+RP
melt_vcf<-readVcf("NA12878_small.melt.vcf")

info(melt_vcf)$support=0.5*(info(melt_vcf)$LP + info(melt_vcf)$RP)
info(melt_vcf)$tot_depth=info(melt_vcf)$LP + info(melt_vcf)$RP


melt.bpr=breakpointRanges(melt_vcf, info_columns=c("support", "tot_depth"))

# assign caller
melt.bpr$Caller="Melt"

# filter
melt.bpr=melt.bpr[which(melt.bpr$FILTER=="PASS")]
