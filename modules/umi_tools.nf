process extract_UMI {
    label 'process_medium'
    tag "$fastq"
    publishDir "${params.output_dir}/fastqs/UMI_extracted",
        mode: 'copy',
        enabled: params.publish_UMI_extracted

    input:
        path fastq
    output:
        path "UMI_$fastq"
    
    script:
    """
    umi_tools extract --bc-pattern=NNNNNN -I $fastq -S UMI_${fastq}
    """
}

process dedup_UMI {
    label 'process_low'
    tag "$bam"

    input:
        path bam
    output:
        path "dedup_$bam", emit: dedup_bams
        path "${bam.simpleName}.log", emit: logs
    
    """
    samtools index $bam
    umi_tools dedup -I $bam -S dedup_$bam > ${bam.simpleName}.log
    """
}