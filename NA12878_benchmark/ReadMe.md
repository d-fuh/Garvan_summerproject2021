---
title: "Benchmarking Structural Variant Callers using WGS data from NA12878"
version: "1.0"
date: "01-14-2022"
---

This task benchmarks the five structural variant callers used (Delly, Manta, Melt, SvABA, and Wham) on Illumina-based WGS data from NA12878.

# Current Working Notes
version 1.0

- Complete the Rmd
- Aim to partition Rmd tasks into individual scripts
- Obtain truth set of NA12878 calls to benchmark complete Section II


# Data processing

Caller outputs (VCF) were processed using [AnnotSV](https://lbgi.fr/AnnotSV/) based on the following criteria:

	# version A
	1. Annotation_mode == "full"
	2. All ACMG classes included ({1:5, NA})
	
	# version B [ONGOING]
	1. Annotation_mode == "split" (Output contain as many annotations lines as the number of genes covered by a given SV)
	2. All ACMG classes included ({1:5, NA})

Regarding the annotation mode, from AnnotSV:

> A typical AnnotSV use would be to first look at the annotation and ranking of each SV as a whole (i.e. “full”) and then focus on the content of that SV. Indeed, there are 2 types of lines produced by AnnotSV (cf the “AnnotSV type” output column):

- An annotation on the “full” length of the SV:
Every SV are reported, even those not covering a gene. This type of annotation gives an estimate of the SV event itself.

- An annotation of the SV “split” by gene:
This type of annotation gives an opportunity to focus on each gene overlapped by the SV. Thus, when a SV spans over several genes, the output will contain as many annotations lines as covered genes (cf example in FAQ). This latter annotation is extremely powerful to shorten the identification of mutation implicated in a specific gene.


# Table of Contents



[TBC]
