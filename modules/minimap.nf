process align_MINIMAP {
    tag "$fastq.simpleName"

    input:
        path ref_fasta
        path fastq
    output:
        path "*.bam"
    
    """
    minimap2 -x map-ont -a --sr -t ${task.cpus} $ref_fasta $fastq -o ${fastq}.sam
    samtools sort -@${task.cpus} -o ${fastq}.bam ${fastq}.sam
    """
}
