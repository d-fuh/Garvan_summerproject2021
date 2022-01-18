---
title: "Benchmarking Structural Variant Callers using WGS data from NA12878"
version: "1.0"
date: "01-17-2022"
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
## 1. SV call location across NA12878
### 1-1. SV location by coverage (small, medium, large)
### 1-2. SV location by callers

## 2. Types of Structural Variants
### 2-1. SV type between callers
### 2-2. SV types vs coverage level

## 3. Size and length of Structural Variants in the call set
### 3-1. SV length by caller
### 3-2. SV length by caller (< 20 kb)
### 3-3. SV length across chromosomes vs caller
### 3-4. SV length vs SV type 

## 4. Number of variants detected by all methods

## 5. ACMG class of variants

## 6. Detected variants affecting CDS
### 6-1. CDS-affecting SVs
### 6-2. Gene counts, exon counts, and frameshift
### 6-3. The "elite list" of variants

## 7. nMDS-ExAC Z scores: A bundled measure of functional constraint
### 7-1. ExACZ ~ Caller
### 7-2. ExACZ ~ ACMG class
### 7-3. ExACZ ~ SV type
### 7-4. ExACZ ~ SV length
### 7-5. ExACZ ~ variant location

## 8. PCA
### 8-1. ExAC Z continued: ACMG score and ExAC Z measure of functional constraint
