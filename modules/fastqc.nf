process fastqc {
    tag "${fastq_files}"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    
    input:
        path fastq_files
    output:
        path "${task.process}/"
    
    """
    mkdir ${task.process}
    fastqc \\
        -t 4 \\
        -o ${task.process}/ \\
        ${fastq_files}
    """
}