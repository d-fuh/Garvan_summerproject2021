The VCF files from the SvABA workflow requires extensive processing to achieve compatibility with other caller outputs from the GATK-SV workflow.

We compiled the following scripts as a pipeline to process SvABA output (VCF) in order to benchmark the results:

### Pipeline (Using MacOS is advised)
#### Input: SvABA *.vcf.gz
#### 1. Use SvABA_header_trimmer.R to remove redundant columns
#### 2. Decompress the output from (1)
#### 3. (Optional) Use SvABA_BND_converter to convert the sv type column (all BND) to proper SV types (INV/DEL/DUP-INS) based on breakpoint orientation
 
The file with sv types converted are stored in ./Converted_csv_output as well as ../VCF

Note that SvABA_further_annotation was not used in our benchmarking workflow and is optional. Source: https://github.com/walaj/svaba/blob/master/R/svaba-annotate.R