process STAR_ALIGN {
    label 'process_max'
    tag "$fastq"
    publishDir "${params.output_dir}/bams",
        mode: 'copy',
        pattern: "*.bam",
        enabled: { params.enable_UMI_treatment ? false : params.publish_mapped }

    input:
        path fastq
        path genome_index

    output:
        path "*.bam", emit: bams
        path "*Log.final.out", emit: logs

    script:
    // remove UMI_ prefix from fastq file name
    outName = fastq.simpleName - "UMI_"
    outName = outName - "trimmed_final_"
    """
    STAR \\
        --runThreadN $task.cpus \\
        --genomeDir $genome_index \\
        --readFilesIn $fastq \\
        --readFilesCommand zcat \\
        --outFileNamePrefix $outName \\
        --outSAMtype BAM SortedByCoordinate \\
        --outReadsUnmapped Fastx \\
        --limitBAMsortRAM 10000000000 \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --outFilterScoreMinOverLread 0.66 \\
        --outFilterMatchNminOverLread 0.66 \\
        --alignEndsType Local \\
        --outSAMattributes Standard \\
        ${params.enable_shared_memory ? '--genomeLoad LoadAndRemove' : ''} 
    """
}