This is where annotated VCF files (in tsv; via AnnotSV) along with original VCF files are stored. 

Scripts to generate analysis-ready data frames are the following:

- maestro_generation_full (for AM=full)
- maestro_generation_split (for AM=split)

The output data frame will be written into the following csv files:

- maestro.full.csv
- maestro.split.csv

..., respectively. These can be read in the main benchmark workflow NA12878_maestro.Rmd. 