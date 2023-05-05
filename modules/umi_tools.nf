process extract_UMI {
    label 'process_medium'
    tag "$fastq"
    publishDir "${params.output_dir}/fastqs/UMI_extracted",
        mode: 'copy',
        enabled: params.publish_UMI_extracted

    input:
        path fastq
    output:
        path "UMI_$fastq"
    
    script:
    def ns = "N" * params.UMI_length
    """
    umi_tools extract --bc-pattern=$ns -I $fastq -S UMI_${fastq}
    """
}

process dedup_UMI {
    label 'process_low'
    tag "$bam"

    input:
        path bam
    output:
        path "dedup_bam/$bam", emit: dedup_bams
        path "${bam.simpleName}.log", emit: logs
    
    """
    samtools index $bam
    mkdir dedup_bam
    umi_tools dedup -I $bam -S dedup_bam/$bam > ${bam.simpleName}.log
    """
}

process cluster_UMI {
    tag "$bam"

    input:
        path bam
    output:
        path "*fastq.gz"

    script:
    def th = params.UMI_threshold
    def window = params.window_size
    """
    # index bam file
    python -c "import pysam; pysam.index('$bam')"
    umiclusterer.py $bam  -j ${task.cpus} -t $th -w $window | gzip > ${bam.simpleName}.fastq.gz 
    """
} 