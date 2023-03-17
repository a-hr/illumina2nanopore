process featureCounts {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.output_dir}",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path dedup_bams
        path raw_bams
        path saf_file
        val prefix
    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs

    """
    # deduplicated bams
    featureCounts -a ${saf_file} -o ${prefix}_dedup_counts.tsv $dedup_bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/Aligned.sortedByCoord.out.bam//g' ${prefix}_dedup_counts.tsv

    # raw bams
    featureCounts -a ${saf_file} -o ${prefix}_raw_counts.tsv $raw_bams -F 'SAF'
    sed -i '1d ; 2 s/Aligned.sortedByCoord.out.bam//g' ${prefix}_raw_counts.tsv
    """
}