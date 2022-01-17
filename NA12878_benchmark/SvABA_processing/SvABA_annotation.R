# Read in svaba vcf file
library(dplyr)

svaba.sv.vcf="~/Documents/GitHub/summerproject2021/NA12878_benchmark/SvABA_processing/NA12878.svaba.vcf"

cols <- colnames(read.table(pipe(paste0('grep -v "##" ', svaba.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
svaba_uniq <- read.table(svaba.sv.vcf, col.names = cols, stringsAsFactors = FALSE)

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

write.csv(svaba_uniq, "svaba.uniq.csv")