process merge_fastqs {
    publishDir "${params.output_dir}/fastqs/basecalled", mode: 'copy', enabled: params.publish_merged

    input:
        val sample
        path "${sample}_??.fastq.gz"
    output:
        path "merged_*.fastq.gz"

    """
    cat ${sample}*.fastq.gz > merged_${sample}.fastq.gz
    """
}