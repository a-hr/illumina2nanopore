process featureCounts_targeted {
    tag "$saf_file"
    label 'process_low'
    publishDir "${params.output_dir}/counts", mode: 'copy', pattern: "*.tsv", enabled: !params.enable_grouping

    input:
        tuple path(saf_file), path(bams)
        val prefix
    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs

    script:
    def library = saf_file.baseName
    """
    featureCounts -T ${task.cpus} -a $saf_file -o ${prefix}_${library}_counts.tsv $bams -F 'SAF'
    sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g' ${prefix}_${library}_counts.tsv
    """

}

process featureCounts_global {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.output_dir}/counts", mode: 'copy', pattern: "*.tsv"

    input:
        path bams
        path gtf_file
        val prefix

    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs
    
    script:

    if (params.fraction_allowed || params.multimapping_allowed) {
        def multimap = params.multimapping_allowed ? ' -M' : ''
        def fraction = params.fraction_allowed ? ' --fraction' : ''

        """
        # create multimapping count table
        featureCounts -T ${task.cpus} ${multimap} ${fraction} -a ${annotations} -o ${prefix}_multi_counts.tsv $bams
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g; 2 s/trimmed_final_//g' ${prefix}_multi_counts.tsv

        # creaet standard count table, allows for contrasts to detect which genes multimapped
        featureCounts -T ${task.cpus} -a ${annotations} -o ${prefix}_counts.tsv $bams
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g; 2 s/trimmed_final_//g' ${prefix}_counts.tsv
        """
    } else {
        """
        featureCounts -T ${task.cpus} -a ${gtf_file} -o ${prefix}_counts.tsv $bams
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g; 2 s/trimmed_final_//g' ${prefix}_counts.tsv
        """
    }
}
