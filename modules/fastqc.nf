process fastqc {
    tag "${fastq_files}"
    
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