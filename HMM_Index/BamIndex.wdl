## source: https://github.com/cancerit/BRASS-wdl/blob/develop/brass-bam-index.wdl

task samtoolsSortTask{
    File inputBAM
    String sortedBAM = basename(inputBAM, ".bam") + ".sorted.bam"
    Int diskSpace = 2*ceil(size(inputBAM, "GB"))
    Int memoryGB
    Int cpu
    Int preemptible
    
    command {
        samtools sort -m 20G -@ 4 -o ${sortedBAM} ${inputBAM} && \
        samtools index ${sortedBAM}
    }
    runtime {
        docker: "erictdawson/svdocker"
        disks: "local-disk ${diskSpace} HDD"
        memory: memoryGB + "GB"
        cpu: cpu
        preemptible: preemptible
    }
    output{
        File outputBAM="${sortedBAM}"
        File outputBAI="${sortedBAM}.bai"
    }
}

workflow samtoolsSort {
    File inputBAM
  
    call samtoolsSortTask{
        input:
            inputBAM=inputBAM
    }
}