## Task 5: Intro to WDL

## Objective

## Inputs:
##   (I1) BAM File (.bam)
##   (I2) BAM Index (.bai)

## Outputs:
##   (O1) File containing the header of the BAM file
##   (O2) File indicating the mapping rate
##   (O3) File indicating the average depth


## v4.0 (final)
## Resolved command syntax bug
## Fixing filename input issue
## Confirmed wf output syntax correct 
## Consider filling in missing runtime parameter (optional)


# START WORKFLOW
version 1.0
workflow WF_HeaderMaprateDepth{
	input{
		File input_bam
		File index_bam
		String filename
		Int preemptible = 3
		Int cpu = 1
		Int memoryGB = 5
		Int diskGB_boot = 5
	}

	call headCaller{
		input:
			input_bam=input_bam,
			index_bam=index_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}

	call mapCaller{
		input:
			input_bam=input_bam,
			index_bam=index_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}

	call depthCaller{
		input:
			input_bam=input_bam,
			index_bam=index_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}
## check the output file format 
	output {
		File header = headCaller.header
		File maprate = mapCaller.maprate
		File meandepth = depthCaller.meandepth
	}
}
# END WORKFLOW
# DEFINE TASKS
# This task calls the header from a BAM file.
task headCaller{
	input{
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}

	command {
		samtools view -H ${input_bam} > ${filename}_header.txt
	}

	output {
		File header = "${filename}_header.txt"
	}

	runtime {
		docker: "staphb/samtools:latest"
    	preemptible: preemptible
    	bootDiskSizeGb: diskGB_boot
    	cpu: cpu
    	memory: memoryGB

	}
}
# This task calls the mapping rate from a BAM file by highlighting the "%" chr.
task mapCaller {
	input {
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}
	
	command {
		samtools flagstat ${input_bam} | grep "%" > ${filename}_maprate.txt
	}

	output {
		File maprate = "${filename}_maprate.txt"
	}

	runtime {
		docker: "staphb/samtools:latest"
    	preemptible: preemptible
    	bootDiskSizeGb: diskGB_boot
    	cpu: cpu
    	memory: memoryGB

	}
}
# This task retrieves the mean depth from a BAM file by means of samtools:coverage.
task depthCaller {
	input {
		File input_bam
		File index_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}
	
	command {
		samtools coverage ${input_bam} | grep -A1 -w "meandepth" > ${filename}_meandepth.txt
	}

	output {
		File meandepth = "${filename}_meandepth.txt"
	}

	runtime {
		docker: "staphb/samtools:latest"
    	preemptible: preemptible
    	bootDiskSizeGb: diskGB_boot
    	cpu: cpu
    	memory: memoryGB

	}
}
