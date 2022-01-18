# Join by ID and plug in the annotated SV_type into old benchmark frames

library(dplyr)

svaba.uniq <- read.csv("SvABA_processing/svaba.uniq.csv", header = TRUE)
benchmark.full <- read.csv("Maestro_scripts/maestro.full.csv")

svaba.uniq$Caller <- "SvABA"
svaba.uniq$SV_length <- 1
#names(svaba.uniq)[names(svaba.uniq) == 'SVTYPE'] <- 'SV_type'


benchmark.full$SV_type[match(svaba.uniq$ID, benchmark.full$ID)] <- svaba.uniq$SV_type
benchmark.full$SV_length[match(svaba.uniq$ID, benchmark.full$ID)] <- svaba.uniq$SV_length
#benchmark.split$SV_type[match(svaba.uniq$ID, benchmark.split$ID)] <- svaba.uniq$SV_type
#benchmark.split$SV_length[match(svaba.uniq$ID, benchmark.split$ID)] <- svaba.uniq$SV_length


# test
benchmark.full %>% filter(Caller == "SvABA") %>%
  ggplot(.) + aes(x=SV_type) + geom_bar()

benchmark.full %>%
  ggplot(.) + aes(x=SV_type, fill=Caller) + geom_bar() +
  theme_bw() + scale_fill_brewer(palette = 2)
