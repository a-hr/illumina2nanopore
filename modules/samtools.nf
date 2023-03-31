process BAM_INDEX {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams",
        mode: 'copy',
        enabled: params.publish_mapped

    input:
        path bam
    output:
        path "*.bam", includeInputs: true
        path "*.bai"
    when:
        params.publish_mapped == true

    """
    samtools index $bam
    """
}