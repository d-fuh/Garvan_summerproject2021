# This workflow converts the BND notation from the original SvABA VCF output
# to proper SV types according to their breakpoint orientation.

# source:

# **NOTE** Use SvABA_header_trimmer.Rmd first 
# to trim the redundant columns out of the raw VCF file

# **NOTE** the output from SvABA_header_trimmer.Rmd
# needs to be decompressed first on MacOS
# regardless of whether the input VCF was decompressed or not

library(dplyr)

svaba.vcf="~/Dropbox/ExRes/ExResVCF/Svaba_backup/NA12878.svaba.vcf"

cols <- colnames(read.table(pipe(paste0('grep -v "##" ', svaba.vcf,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))

svaba_uniq <- read.table(svaba.vcf, col.names = cols, stringsAsFactors = FALSE)

# Assign SV type
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

# write the converted object to csv
# the sv type info can be extracted as a single column
# then pasted to VCF/GR objects.
write.csv(svaba_uniq, "NA12878_svaba_converted.csv")

# can also write into VCF files directly?
