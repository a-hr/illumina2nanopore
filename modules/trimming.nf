process adapter_trim {
    // trims the adapters, filters by length and demultiplexes by strand (forward/reverse)
    tag "${fastq_file}"
    publishDir "${params.output_dir}/fastqs/trimmed_fastqs", pattern: "trimmed_*.fastq.gz", mode: 'copy'
    
    input:
        tuple val(five_prime_fw), val(three_prime_fw)
        tuple val(five_prime_rv), val(three_prime_rv)
        path fastq_file
    output:
        path "trimmed_*.fastq.gz", emit: trimmed_fastqs
        path "*_demultiplex.log", emit: log
    
    """
    cutadapt \\
        -j 0 \\
        -a forward="${five_prime_fw}...${three_prime_fw}" \\
        -a reverse="${five_prime_rv}...${three_prime_rv}" \\
        -e 0.2 \\
        -m 50 \\
        -M 300 \\
        -o trimmed_{name}.fastq.gz \\
        ${fastq_file} > ${fastq_file.simpleName}_demultiplex.log
    """
}