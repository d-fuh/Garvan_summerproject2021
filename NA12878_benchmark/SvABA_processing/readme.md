This is the directory for raw SvABA VCF file processing workflow.

# Workflow Overview

SvABA VCF generation -> svaba.sv.vcf -> SvABA_output_processing.rmd 

-> svaba.new.vcf -> SvABA_annotation.R -> svaba.uniq (R data frame/csv)

-> (optional) SvABA_further_annotation [INCOMPLETE]

-> merge_svaba_to_benchmark.R -> benchmark.full/split (R data frame/csv)

-> NA12878_benchmark/maestro2.0.Rmd (main benchmark workflow)
