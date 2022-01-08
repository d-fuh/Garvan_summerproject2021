---
title: "Benchmarking Structural Variant Callers using WGS data NA12878"
date: "01-08-2022"
---
***Version 1.0***

This task benchmarks the five structural variant callers used (Delly, Manta, Melt, SvABA, and Wham) on the WGS data NA12878.

# Data processing

Caller outputs (.vcf.gz) were processed using [AnnotSV](https://lbgi.fr/AnnotSV/) based on the following criteria:

	1. Annotation_mode == "split" (Output contain as many annotations lines as genes covered by SV)
	2. Rank = {4,5} (Excluding NA, and any calls that were ACMG-classified as non-pathogenic)

For future editions: Investigate how to subset using "the number of supporting reads".

Regarding the annotation mode, from AnnotSV:

> A typical AnnotSV use would be to first look at the annotation and ranking of each SV as a whole (i.e. “full”) and then focus on the content of that SV. Indeed, there are 2 types of lines produced by AnnotSV (cf the “AnnotSV type” output column):

- An annotation on the “full” length of the SV:
Every SV are reported, even those not covering a gene. This type of annotation gives an estimate of the SV event itself.

- An annotation of the SV “split” by gene:
This type of annotation gives an opportunity to focus on each gene overlapped by the SV. Thus, when a SV spans over several genes, the output will contain as many annotations lines as covered genes (cf example in FAQ). This latter annotation is extremely powerful to shorten the identification of mutation implicated in a specific gene.
