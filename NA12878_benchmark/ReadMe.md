---
title: "Benchmarking Structural Variant Callers using WGS data NA12878"
date: "01-08-2022"
---
***Version 1.0***

This task benchmarks the five structural variant callers used (Delly, Manta, Melt, SvABA, and Wham) on the WGS data NA12878.

# Data processing

Caller outputs (.vcf.gz) were processed using [AnnotSV](https://lbgi.fr/AnnotSV/) based on the following criteria:

	1. Annotation_mode == "split" (Only
	2. Rank = {4,5} (Excluding NA, and any calls that were ACMG-classified as non-pathogenic)

For future editions: Investigate how to subset using "the number of supporting reads".

