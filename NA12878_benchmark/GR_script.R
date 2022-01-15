# Generate GRange object from VCFs

library(vcfR)
library(GenomicRanges)

path.expand('~')
# omittable
delly.vcf <- read.vcfR("./GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/VCF/NA12878.delly.vcf")

print(delly.vcf@fix)

# can we generate gr from annot.tsv?
#benchmark.full <- read.csv("./GitRepo/summerproject2021/NA12878_benchmark/Maestro_scripts/maestro.full.csv")

large.tst <- large.full %>% 
  subset(., select = c("SV_chrom", "SV_start", "SV_end", "QUAL", "Caller", "SV_type"))

# generate gr from delly.vcf as an example
## delly
filter(large.tst, Caller == "Delly") -> delly.t

delly.gr <- GRanges(
  seqnames = Rle(delly.t$SV_chrom, delly.t$rownames),
  ranges = IRanges(delly.t$SV_start, end = delly.t$SV_end, names = delly.t$rownames),
)

## manta
filter(large.tst, Caller == "Manta") -> manta.t

manta.gr <- GRanges(
  seqnames = Rle(manta.t$SV_chrom, manta.t$rownames),
  ranges = IRanges(manta.t$SV_start, end = manta.t$SV_end, names = manta.t$rownames),
)

## svaba
filter(large.tst, Caller == "SvABA") -> svaba.t

svaba.gr <- GRanges(
  seqnames = Rle(svaba.t$SV_chrom, svaba.t$rownames),
  ranges = IRanges(svaba.t$SV_start, end = svaba.t$SV_end, names = svaba.t$rownames),
)

## melt
filter(large.tst, Caller == "Melt") -> melt.t

melt.gr <- GRanges(
  seqnames = Rle(melt.t$SV_chrom, melt.t$rownames),
  ranges = IRanges(melt.t$SV_start, end = melt.t$SV_end, names = melt.t$rownames),
)

## wham
filter(large.tst, Caller == "Wham") -> wham.t

wham.gr <- GRanges(
  seqnames = Rle(wham.t$SV_chrom, wham.t$rownames),
  ranges = IRanges(wham.t$SV_start, end = wham.t$SV_end, names = wham.t$rownames),
)

# find overlaps
findOverlaps(delly.gr, manta.gr)

