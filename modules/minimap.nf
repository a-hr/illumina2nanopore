process MINIMAP_ALIGN {
    input:
        path ref_fasta
        path fastq
    output:
        path "*.bam"
    
    """
    minimap2 -x map-ont -a --sr -t ${task.cpus} $ref_fasta $fastq | samtools sort @ ${task.cpus}-o ${fastq.simpleName}.bam  
    """
}
