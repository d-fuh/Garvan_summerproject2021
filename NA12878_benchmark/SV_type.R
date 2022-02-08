# BND from csv to svaba_svgr by 
small.svaba.uniq=read.csv("~/Dropbox/NA12878/small.svaba.uniq.csv")
med.svaba.uniq=read.csv("~/Dropbox/NA12878/med.svaba.uniq.csv")
ft.svaba.uniq=read.csv("~/Dropbox/NA12878/small_50pc_converted.csv")
large.svaba.uniq=read.csv("~/Dropbox/NA12878/svaba.uniq.csv")

# svgr#elementMetadata$sourceID == csv$ID
# svgr@elementMetadata$svtype == csv$SV_type

small_svaba_svgr@elementMetadata$svtype[which(small_svaba_svgr@elementMetadata$sourceId == small.svaba.uniq$ID)]=small.svaba.uniq$SV_type

med_svaba_svgr@elementMetadata$svtype[which(med_svaba_svgr@elementMetadata$sourceId == med.svaba.uniq$ID)]=med.svaba.uniq$SV_type

ft_svaba_svgr@elementMetadata$svtype[which(ft_svaba_svgr@elementMetadata$sourceId == ft.svaba.uniq$ID)]=ft.svaba.uniq$SV_type

svaba_svgr@elementMetadata$svtype[which(svaba_svgr@elementMetadata$sourceId == large.svaba.uniq$ID)]=large.svaba.uniq$SV_type

# sv type
dflarge=large_svgr@elementMetadata[, c("svtype", "Caller")]
dfft=ft_svgr@elementMetadata[, c("svtype", "Caller")]
dfmed=med_svgr@elementMetadata[, c("svtype", "Caller")]
dfsmall=small_svgr@elementMetadata[, c("svtype", "Caller")]

# cov
dflarge$Cov="35x"
dfft$Cov="18x"
dfmed$Cov="6x"
dfsmall$Cov="1x"

#
dflarge=as.data.frame(dflarge)
dfft=as.data.frame(dfft)
dfmed=as.data.frame(dfmed)
dfsmall=as.data.frame(dfsmall)

# join
plyr::join(dflarge, dfft, type="full") %>% plyr::join(dfmed, type="full") %>% plyr::join(dfsmall, type="full") -> svt

# plot
ggplot(svt) + aes(x=svtype, fill=Caller) + geom_bar(position="dodge") + facet_wrap(~Cov) + theme_bw()
ggsave("SV_type_by_coverage.png")

# reorder panel header level for PR plot
svt$Cov <- factor(svt$Cov, levels = c("35x", "18x", "6x", "1x"))
