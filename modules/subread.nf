process featureCounts {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.output_dir}",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path bams
        path annotations
        val prefix
    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs

    script:
    def multimap = params.multimapping_allowed ? ' -M' : ''
    def fraction = params.fraction_allowed ? ' --fraction' : ''

    if (params.enable_isoform_counting) {
        """
        featureCounts -T ${task.cpus} ${multimap} ${fraction} -a ${annotations} -o ${prefix}_counts.tsv $bams -F 'SAF'
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g' ${prefix}_counts.tsv
        """
    }
    else {
        """
        featureCounts -T ${task.cpus} ${multimap} ${fraction} -a ${annotations} -o ${prefix}_counts.tsv $bams
        """
    }
}