process MINIMAP_ALIGN {
    tag "$fastq.simpleName"
    publishDir "${params.output_dir}/bams",
        mode: 'copy',
        pattern: "*.bam",
        enabled: { params.enable_UMI_treatment ? false : params.publish_mapped }

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
