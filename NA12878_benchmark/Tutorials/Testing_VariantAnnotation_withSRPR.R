## Testing David's functions

setwd("~/GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF")

library(StructuralVariantAnnotation)

#### Hack the bed file to convert to VCF: save as a txt file and append the correct header in terminal
# personalis <- read.table("Spiral_Genetics_insertions_hg38.bed", header=F)
# personalis$ID=paste0("DUP", seq(1:nrow(personalis)))
# personalis$REF="N"
# personalis$ALT="<DUP>"
# personalis$QUAL="."
# personalis$FILTER="PASS"
# personalis$INFO=paste0("END=",personalis$V3,";SVLEN=", personalis$V3-personalis$V2, 
#                       ";SVTYPE=DUP")
# colnames(personalis)[1:3]=c("CHROM", "POS", "END")
# write.table(personalis[ ,-3], "~/Downloads/test_spiral.txt", sep="\t", row.names = F, quote = F,
#             col.names = F)
# create GR object

## Load Example truth set: personalis

truth_personalis <- readVcf("test_personalis_hg38.vcf")
truth_personalis2=breakpointRanges(truth_personalis)

## Example: Load data from delly, manta and svaba

delly_vcf<-readVcf("NA12878_small.delly.vcf.gz")

head(delly_vcf)
## Note that there are 2 main components of this object:
# 1. an "info" data.frame which contains the chr location + information in "INFO" column
# 2. a "geno" list which contains information in the "GT column - this is also the NA12878 genotype information"

## To Get the information on Paired Reads and Spanning Reads
## Note that the header of Delly says:
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">									
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">	
## This info is stored in INFO, so we need to pull it out using that command
head(info(delly_vcf))
# you can pull this info into the resulting breakpointRanges object using "info_columns"
delly_vcf2=breakpointRanges(delly_vcf, info_columns=c("SR", "PE"))
head(delly_vcf2)
# we can now append this data to delly_vcf2
delly_vcf2$PE=PE
delly_vcf2$SR=SR

## Assign a column for caller name
delly_vcf2$caller="delly"
## Need to filter out the variants that pass
delly_vcf3=delly_vcf2[which(delly_vcf2$FILTER=="PASS")]

## Do the same for svaba
## 1. Need to delete the columns with no header
## 2. Need to run the svaba_annotate_script.R
## Was there a file to change the annotations for svaba from BND to DUP/INS etc
svaba_vcf<-readVcf("NA12878_svaba_troubleshoot/test.vcf")
svaba_vcf2=breakpointRanges(svaba_vcf)
svaba_vcf2$caller="svaba"
## See if you can insert the svaba_annotate_script here to reassign the svaba_vcf2$svtype.
## You can load the file separately and save SVTYPE as a new variable or to file

## Note that the header of Svaba says:
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of discordant-supported reads for this variant">												
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of spanning reads for this variants">													
## This information is in FORMAT and NOT THE INFO column, 
## We need to extract this from the GT or genotype info, and save it to the INFO data frame
head(geno(svaba_vcf))
SR=geno(svaba_vcf)$SR
head(SR)
## If the 3 columns are not deleted, there should be a matrix with 4 columns - we only want the final column
info(svaba_vcf)$SR=SR[ ,4]
head(info(svaba_vcf))
# Now, do the same thing for the paired or discordant reads!!
PE=???

## find overlaps between delly and svaba and map these to the delly file:
## These are saved in a column called svabamatch
## For this you may want to find the file that has the largest number of variants and test this
## Or you can use a combined list 
delly_vcf3$svabamatch=countBreakpointOverlaps(delly_vcf3, svaba_vcf2, maxgap=100, sizemargin = 0.25, restrictMarginToSizeMultiple = 0.5, countOnlyBest=T)
head(delly_vcf3)
table(delly_vcf3$svabamatch)

## find overlaps between svaba and truth
truth_personalis2$svabamatch=countBreakpointOverlaps(truth_personalis2, svaba_vcf2, maxgap=100, sizemargin = 0.25, restrictMarginToSizeMultiple = 0.5, countOnlyBest=T)
head(truth_personalis2)
table(truth_personalis2$svabamatch)

## we can also do the reverse and annotate your delly_vcf2 with which variants overlap with the truth set.
## This way
delly_vcf3$truthmatch=countBreakpointOverlaps(delly_vcf3, truth_personalis2, maxgap=100, sizemargin = 0.25, restrictMarginToSizeMultiple = 0.5, countOnlyBest=T)
head(delly_vcf3)
# now you can try ti compare SR and PE for the samples which are matched to the truthset
summary(delly_vcf3$SR[which(delly_vcf3$truthmatch==1)])