# This workflow converts the BND notation from SvABA output
# to appropriate SV types in $SV_type

# Use SvABA_header_trimmer.Rmd first 
# to trim the redundant columns
# in the original vcf output

# Note that the output from SvABA_header_trimmer.Rmd
# needs to be decompressed first
# regardless of whether the input vcf was decompressed
# for that workflow

library(dplyr)

svaba.vcf="~/Dropbox/ExRes/ExResVCF/Svaba_backup/ER019_SAR1N.svaba.vcf"

cols <- colnames(read.table(pipe(paste0('grep -v "##" ', svaba.vcf,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))

svaba_uniq <- read.table(svaba.vcf, col.names = cols, stringsAsFactors = FALSE)

get_sv_type <- function(x){
  # Find mate pair
  root <- gsub(":[12]", "", x)
  mate1 <- paste0(root, ":1")
  mate2 <- paste0(root, ":2")
  alt1 <- svaba_uniq %>% filter(ID == mate1) %>% .$ALT
  alt2 <- svaba_uniq %>% filter(ID == mate2) %>% .$ALT
  # Determine sv type based on breakpoint orientation
  if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))){
    sv_type <- "INV"
    
  } else if (grepl("[a-zA-Z]\\[", alt1) & grepl("^\\]", alt2)){
    sv_type <- "DEL"
    
  } else if (grepl("^\\]", alt1) & grepl("[a-zA-Z]\\[", alt2)){
    sv_type <- "DUP/INS"
    
  } else{
    sv_type <- "UNKNOWN"
  }
  return(sv_type)
}

svaba_uniq$SV_type <- sapply(svaba_uniq$ID, get_sv_type)

# write converted file to csv
# because write.vcf doesn't work
# then extract in your following workflow
# the $SV_type column to your new vcf
write.csv(svaba_uniq, "ER019_SAR1N.svtype.csv")
