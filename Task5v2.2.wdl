## Task 5: Intro to WDL

## Objective

## Inputs:
##   (I1) BAM File (.bam)
##   (I2) BAM Index (.bai)

## Outputs:
##   (O1) File containing the header of the BAM file
##   (O2) File indicating the mapping rate
##   (O3) File indicating the average depth


## Workflow construction
## v2.2 The no-mount version
version 1.0
# DEFINE WORKFLOW
workflow WF_HeaderMaprateDepth{

	input{
		File input_bam = input_bam
		String filename = filename
		Int preemptible = 3
		Int cpu = 1
		Int memoryGB = 5
		Int diskGB_boot = 5
	}

	call headCaller{
		input:
			input_bam=input_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}

	call mapCaller{
		input:
			input_bam=input_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}

	call depthCaller{
		input:
			input_bam=input_bam,
			filename=filename,
			memoryGB=memoryGB,
			diskGB_boot=diskGB_boot,
			preemptible=preemptible,
			cpu=cpu
	}
}

# DEFINE TASKS
task headCaller{
	input{
		File input_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}

	command <<<
	# grep -w "ID" for Unix

		samtools view -H ${input_bam} | Select-String -Pattern "ID" > ${filename}_header.txt
	>>>

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

task mapCaller {
	input {
		File input_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}
	
	command <<<
	# grep "%" for Linux

		samtools flagstat ${input_bam} | Select-String -Pattern "%"  > ${filename}_maprate.txt
	>>>

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

task depthCaller {
	input {
		File input_bam
		String filename
		Int preemptible
		Int diskGB_boot
		Int cpu
		Int memoryGB
	}
	
	command <<<
	# grep -A1 -w "meandepth" for Unix

		samtools coverage ${input_bam} | Select-String -Pattern "meandepth" -Context 0,1 > ${filename}_meandepth.txt
	>>>

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

## Remarks for future version
## 1. Unsure if .bai needed
## 2. Unsure if command <<< >>> valid syntax
## 3. Unsure if system configuration specified correctly & how it works