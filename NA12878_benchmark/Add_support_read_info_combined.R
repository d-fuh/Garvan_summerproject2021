# support read info extraction script
# for NA12878 small subset
library(StructuralVariantAnnotation)
setwd("~/GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF")


# Below are the stats we want to extract

# |Caller|Original stats used|"Combined PE+SR score" | Note                 |
# |:-----|:------------------|:---------------------:|:---------------------|
# |Delly |SR, PE, RC (X)     | DR                    |                      |
# |SvABA |SR, DR, AD         | AD                    | DP=total depth       |
# |Manta |SR, PR, BND_DEPTH  | SR+ 1/2PR?            | SR lists both the alt & ref so you need to split it  |
# |Wham  |SR, SP             | A                     |                      |
# |Melt  |LP, RP             | LP + RP?              |                      |



## Delly: SR, PE, DR
delly_vcf<-readVcf("NA12878_small.delly.vcf")

# Both SR, PR are already in INFO
# DR is in geno so import to INFO
info(delly_vcf)$DR=geno(delly_vcf)$DR
##can check
##rownames(geno(delly_vcf)$DR)==rownames(info(delly_vcf))

# create svgr object
delly.bpr=breakpointRanges(delly_vcf, info_columns=c("SR", "PE", "DR"))

# Assign a column for caller name
delly.bpr$Caller="Delly"

# Need to filter out the variants that pass
delly.bpr=delly.bpr[which(delly.bpr$FILTER=="PASS")]




## Svaba: SR, DR, AD
svaba_vcf<-readVcf("NA12878_small.svaba.vcf")

# original svaba vcf has $SVTYPE == "BND" for all calls
# used a script to convert them into normal SV_type
# and written into svaba_converted.csv files
small.conv.svaba <- read.csv("small_svaba_converted.csv") %>% subset(select = -X)
info(svaba_vcf)$SVTYPE=small.conv.svaba$SV_type ## double check the ALT column matches
    #svaba_vcf@fixed$ALT == small.conv.svaba$ALT

# import support read info
info(svaba_vcf)$AD=geno(svaba_vcf)$AD
info(svaba_vcf)$DR=geno(svaba_vcf)$DR
info(svaba_vcf)$SR=geno(svaba_vcf)$SR

# check rowname ID matches
  #rownames(AD) == rownames(info(svaba_vcf))

# create svgr object
svaba.bpr<-breakpointRanges(svaba_vcf, info_columns=c("SR", "DR", "AD"))

# Assign a column for caller name
svaba.bpr$Caller="SvABA"

# Need to filter out the variants that pass
svaba.bpr=svaba.bpr[which(svaba.bpr$FILTER=="PASS")]




## Manta: PR, SR, BND_DEPTH, **NEW** SR+0.5*PR
manta_vcf<-readVcf("NA12878_small.manta.vcf")
    #head(manta_vcf)

# find PR, SR, BND_DEPTH
  #colnames(info(manta_vcf))
  #colnames(geno(manta_vcf))

# PR, SR are in GT; BND_DEPTH is in INFO
# Import PR, SR into INFO
  info(manta_vcf)$PR=geno(manta_vcf)$PR ## does PR=PE?
  info(manta_vcf)$SR=geno(manta_vcf)$SR

# check rowname ID matches
  #rownames(geno(manta_vcf)$PR) == rownames(info(manta_vcf))
  #rownames(geno(manta_vcf)$SR) == rownames(info(manta_vcf))

# create svgr object
  # don't forget BND_DEPTH in INFO
manta.bpr=breakpointRanges(manta_vcf, info_columns=c("PR", "SR", "BND_DEPTH"))
  #head(manta.bpr) # BND_DEPTH has lots of NAs as-is

# want to mutate X=SR+0.5*PR
# but dplyr doesn't work with S4 objects
# use base function instead
manta.bpr$X = manta.bpr$SR + 0.5*manta.bpr$PR
  ## Error in 0.5 * manta.bpr$PR : non-numeric argument to binary operator
  ## both variables are in "list" format

manta.bpr$X = unlist(manta.bpr$SR) + 0.5*unlist(manta.bpr$PR)
  ##Error in `[[<-`(`*tmp*`, name, value = c(0, 4.5, 0, 3.5, 19.5, 9, 8.5,  : 
  ##676 elements in value to replace 338 elements
  ##In addition: Warning message:
  ##In unlist(manta.bpr$SR) + 0.5 * unlist(manta.bpr$PR) :
  ##longer object length is not a multiple of shorter object length
# Both SR/PR were recorded as 2-tuples: REF - ALT

# Start from the vcf level instead (of svgr level)
# Assign vector to hold this variable

X=as.numeric(unlist(info(manta_vcf)$SR)) + 0.5*as.numeric(unlist(info(manta_vcf)$PR))
  ## works, but the resulting vector is the two uncoupled list added together
  ## and the order is unknown??
  ## to illustrate this:
A=as.numeric(unlist(info(manta_vcf)$SR))
B=0.5*as.numeric(unlist(info(manta_vcf)$PR))

A+B == X

# assign caller
manta.bpr$Caller="Manta"

# filter
manta.bpr=manta.bpr[which(manta.bpr$FILTER=="PASS")]




## Wham: SP, SR
wham_vcf<-readVcf("NA12878_small.wham.vcf")
  #head(wham_vcf)

# find SP, SR
  #colnames(info(wham_vcf)) ## only SR
  #colnames(geno(wham_vcf)) ## gives NULL !? but they're there...
    #head(geno(wham_vcf)) ## 3 columns INCLUDING SP

# extract SP from GT to INTO
info(wham_vcf)$SP=geno(wham_vcf)$SP

# check rowname ID matches
  #rownames(geno(wham_vcf)$SP) == rownames(info(wham_vcf))

# create svgr object
# don't forget SR in INFO
wham.bpr=breakpointRanges(wham_vcf, info_columns=c("SP", "SR"))
  #head(wham.bpr) ## BND_DEPTH has lots of NAs originally

# assign caller
wham.bpr$Caller="Wham"

# filter
wham.bpr=wham.bpr[which(wham.bpr$FILTER=="PASS")]




## Melt: LP, RP
melt_vcf<-readVcf("NA12878_small.melt.vcf")
    #head(melt_vcf)

# find LP, RP
    #colnames(info(melt_vcf)) ## LP, RP

# since LP RP are in INFO already, we can directly import to svgr by
melt.bpr=breakpointRanges(melt_vcf, info_columns=c("LP", "RP"))
    #head(melt.bpr) # BND_DEPTH has lots of NAs originally

# assign caller
melt.bpr$Caller="Melt"

# filter
melt.bpr=melt.bpr[which(melt.bpr$FILTER=="PASS")]
