library(vcfR)
library(GenomicRanges)
library(StructuralVariantAnnotation)

### Manta: Using vcfR to look at your file
MantaIn=read.vcfR("NA12878_small.manta.vcf")
MantaIn=MantaIn
PR=extract.gt(MantaIn, "PR")
SR=extract.gt(MantaIn, "SR")
BND_DEPTH=extract.info(MantaIn, "BND_DEPTH")
GT=extract.gt(MantaIn, "GT")
tempdf=data.frame(PR=PR, SR=SR, BND_DEPTH=BND_DEPTH, GT=GT)
colnames(tempdf)=c("PR", "SR", "BND_DEPTH", "GT")
tempdf[c(1:2, 27), ]
## Note: The first row has ALT 7 PR, 2 SR supporting the depth. 
### So BND_DEPTH might not be a good measure of allelic depth here.
## So you could probably use:
## Depth = total of all PR + SR (e.g. 16 in the 2nd row for MantaBND:2978:0:1:0:0:0:0)
## Depth = 20+7+39+13 [MantaDUP:TANDEM:1292:0:1:0:0:0]
## Reads supporting Variant = PR[2] +SR[2] (ie. 7 in the first row for MantaBND:2978:0:1:0:0:0:0)
## Reads supporting variant = 7 + 13 (MantaDUP:TANDEM:1292:0:1:0:0:0)



## Now to import using readVcf
## Manta: PR, SR, BND_DEPTH, **NEW** SR+0.5*PR

manta_vcf<-readVcf("NA12878_small.manta.vcf")
PR=geno(manta_vcf)$PR
PR[1]
## Note this entire thing is a list, You will need to extract every first item in this list.
PR_ref=sapply(PR, function(x) x[1])
PR_alt=sapply(PR, function(x) x[2])

## Compare these numbers to the original vcfR file:
tempdf$PR_ref=PR_ref
tempdf$PR_alt=PR_alt



### Melt read in?
melt_vcf<-readVcf("NA12878_small.melt.vcf")
supportingReads=info(melt_vcf)$LP+ info(melt_vcf)$RP
