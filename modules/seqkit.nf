process stats {
    tag "${fastq_files} stats"
    publishDir "${params.output_dir}/seqkit", mode: 'copy'
    
    input:
        path fastq_files
        val name
    output:
        path "*_stats.tsv"
    
    """
    seqkit stats -a -T *fastq.gz > ${name}_${fastq_file[0].simpleName}_demultiplex_stats.tsv
    """
}

process reverse_complement {
    tag "${fastq_file.simpleName}"
    publishDir "${params.output_dir}/fastqs/complemented", mode: 'copy', enabled: params.publish_complemented
    
    input:
        path fastq_file
    output:
        path "*.fastq.gz", includeInputs: true
    
    script:
    f1 = "${fastq_file.simpleName}_R1.fastq.gz"
    f2 = "${fastq_file.simpleName}_R2.fastq.gz"
    """
    mv ${fastq_file} ${f1}
    seqkit seq -r -p ${f1} -o ${f2}
    """
}