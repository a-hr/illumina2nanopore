process merge_fastqs {
    publishDir "${params.output_dir}/fastqs/basecalled",
        mode: 'copy',
        enabled: params.publish_merged

    input:
        val sample
        path "${sample}_??.fastq.gz"
    output:
        path "merged_*.fastq.gz"

    """
    cat ${sample}*.fastq.gz > merged_${sample}.fastq.gz
    """
}

process merge_forward_reverse {
    publishDir "${params.output_dir}/fastqs/complemented",
        mode: 'copy',
        enabled: params.publish_merged_complemented

    input:
        path fastqs
    output:
        path "*.fastq.gz"

    """
    cat *.fastq.gz > final_forward.fastq.gz
    """
}