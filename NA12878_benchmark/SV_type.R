# This script imports SV type data from the BND converter (See SvABA_BND_annotation.R)
# to svaba_svgr objects
# and then analyse the type of SVs detected in each caller 
# using NA12878 at different coverage levels

setwd("~/GitRepo/summerproject2021/NA12878_benchmark/SvABA_processing")

# BND from csv to svaba_svgr by 
med.svaba.uniq=read.csv("med_svaba_converted.csv")
ft.svaba.uniq=read.csv("50pc_converted.csv")
large.svaba.uniq=read.csv("svaba_converted.csv")

# svgr#elementMetadata$sourceID ~ csv$ID
# svgr@elementMetadata$svtype ~ csv$SV_type

med_svaba_svgr@elementMetadata$svtype[which(med_svaba_svgr@elementMetadata$sourceId == med.svaba.uniq$ID)]=med.svaba.uniq$SV_type

ft_svaba_svgr@elementMetadata$svtype[which(ft_svaba_svgr@elementMetadata$sourceId == ft.svaba.uniq$ID)]=ft.svaba.uniq$SV_type

svaba_svgr@elementMetadata$svtype[which(svaba_svgr@elementMetadata$sourceId == large.svaba.uniq$ID)]=large.svaba.uniq$SV_type

# sv type
dflarge=large_svgr@elementMetadata[, c("svtype", "Caller")]
dfft=ft_svgr@elementMetadata[, c("svtype", "Caller")]
dfmed=med_svgr@elementMetadata[, c("svtype", "Caller")]

# cov
dflarge$Cov="35X"
dfft$Cov="18X"
dfmed$Cov="6X"


dflarge=as.data.frame(dflarge)
dfft=as.data.frame(dfft)
dfmed=as.data.frame(dfmed)

# join
plyr::join(dflarge, dfft, type="full") %>% plyr::join(dfmed, type="full") -> svt

# plot
ggplot(svt) + aes(x=svtype, fill=Caller) + geom_bar(position="dodge") + facet_wrap(~Cov) + theme_bw()
ggsave("~/GitRepo/summerproject2021/NA12878_benchmark/Figures/SV_type_by_coverage.png")

# reorder panel header level for PR plot
svt$Cov <- factor(svt$Cov, levels = c("35X", "18X", "6X"))
