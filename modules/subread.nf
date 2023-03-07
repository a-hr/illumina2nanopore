process featureCounts {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.output_dir}",
        mode: 'copy',
        pattern: "*.tsv",
        stageAs: {filename -> "${prefix}_${filename}"}

    input:
        path dedup_bams
        path raw_bams
        path saf_file
    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs

    script:
    prefix = dedup_bams[0].name.tokenize('_')[0]
    """
    # deduplicated bams
    featureCounts -a ${saf_file} -o dedup_counts.tsv $dedup_bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/Aligned.sortedByCoord.out.bam//g' dedup_counts.tsv

    # raw bams
    featureCounts -a ${saf_file} -o raw_counts.tsv $raw_bams -F 'SAF'
    sed -i '1d ; 2 s/Aligned.sortedByCoord.out.bam//g' raw_counts.tsv
    """
}