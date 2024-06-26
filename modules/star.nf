process align_STAR {
    label 'process_max'
    tag "$fastq"

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
        --outSAMattributes NH HI AS nM NM MD \\
        ${params.enable_shared_memory ? '--genomeLoad LoadAndRemove' : ''} 
    """
}