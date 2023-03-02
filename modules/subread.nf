process featureCounts {
    label 'process_low'
    publishDir "${params.output_dir}",
        mode: 'copy',
        pattern: "*.txt"

    input:
        path bams
        path saf_file
    output:
        path "counts.txt", emit: counts
        path "*.summary", emit: logs

    """
    featureCounts -a ${saf_file} -o counts.txt $bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/_UMIAligned.sortedByCoord.out.bam//g' counts.txt
    """
}