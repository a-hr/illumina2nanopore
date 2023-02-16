process demultiplex_stats {
    tag "${fastq_file}"
    publishDir "${params.output_dir}/seqkit", mode: 'copy'
    
    input:
        path fastq_file
    output:
        path "*_stats.tsv"
    
    """
    seqkit stats -a -T *fastq.gz > demultiplex_stats.tsv
    """
}

process reverse_complement {
    tag "${fastq_file}"
    publishDir "${params.output_dir}/fastqs/complemented", mode: 'copy'
    
    input:
        path fastq_file
    output:
        path "*.fastq.gz", includeInputs: true
    
    """
    mv ${fastq_file} ${fastq_file.simpleName}_R1.fastq.gz
    seqkit seq -r -p ${fastq_file} | gzip -c  > ${fastq_file.simpleName}_R2.fastq.gz
    """
}