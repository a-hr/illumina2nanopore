process featureCounts {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.output_dir}/counts", mode: 'copy', pattern: "*.tsv", enabled: !params.enable_grouping

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
        def saf = params.saf_files[bams[0].split("_")[0]][0]
        """
        featureCounts -T ${task.cpus} $multimap $fraction -a $saf -o ${prefix}_counts.tsv $bams -F 'SAF'
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g' ${prefix}_counts.tsv
        """
    }

    if (params.fraction_allowed || params.multimapping_allowed) {
        """
        # create multimapping count table
        featureCounts -T ${task.cpus} ${multimap} ${fraction} -a ${annotations} -o ${prefix}_multi_counts.tsv $bams
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g; 2 s/trimmed_final_//g' ${prefix}_multi_counts.tsv

        # creaet standard count table, allows for contrasts to detect which genes multimapped
        featureCounts -T ${task.cpus} -a ${annotations} -o ${prefix}_counts.tsv $bams
        sed -i '1d ; 2 s/${prefix}_//g ; 2 s/Aligned.sortedByCoord.out.bam//g; 2 s/trimmed_final_//g' ${prefix}_counts.tsv
        """
    } 

}