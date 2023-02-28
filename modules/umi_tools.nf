process extract_UMI {
    label 'process_medium'
    tag "$fastq"
    publishDir "${params.output_dir}/fastqs/UMI_extracted",
        mode: 'copy',
        enabled: params.publish_UMI_extracted

    input:
        path fastq
    output:
        path "*_UMI.*"
    
    """
    umi.py -m extract -I $fastq
    """
}

process dedup_UMI {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams",
        mode: 'copy',
        pattern: "dedup_*",
        enabled: params.publish_mapped

    input:
        path bam
    output:
        path "dedup_*", emit: dedup_bams
        path "${logname}.log", emit: logs
    
    script:
    logname = bam.simpleName
    """
    samtools index $bam
    umi_tools dedup -I $bam -S dedup_$bam > ${logname}.log
    """
}