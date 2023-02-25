process adapter_trim {
    // trims the adapters, filters by length and demultiplexes by strand (forward/reverse)
    tag "${fastq_file}"
    publishDir "${params.output_dir}/fastqs/trimmed_fastqs", pattern: "trimmed_*.fastq.gz", mode: 'copy', enabled: params.publish_trimmed
    
    input:
        tuple val(five_prime_fw), val(three_prime_fw)
        tuple val(five_prime_rv), val(three_prime_rv)
        path fastq_file
    output:
        path "trimmed_*.fastq.gz", emit: trimmed_fastqs
        path "*_demultiplex.log", emit: log
    
    // minimum expected seq length after trimming is 145bp
    // maximum expected seq length after trimming is 152bp
    // leave a 5% margin of error in length
    // 152 * 1.05 = 159.6
    // 145 * 0.95 = 137.25

    """
    cutadapt \\
        -j 0 \\
        -a forward="${five_prime_fw}...${three_prime_fw}" \\
        -a reverse="${five_prime_rv}...${three_prime_rv}" \\
        -e 0.2 \\
        -m 137 \\
        -M 160 \\
        -o trimmed_{name}.fastq.gz \\
        ${fastq_file} > ${fastq_file.simpleName}_demultiplex.log
    """
}

process demultiplex_bc {
    tag { prefix }
    publishDir "${params.output_dir}/demultiplexed", mode: 'copy', enabled: params.publish_demultiplexed

    input:
        path fastqs
        path barcode_file

    output:
        path "*.gz", emit: fastqs
        path "*.log", emit: log

    script:
    name = fastqs[0].simpleName
    prefix = name.endsWith('_R1') ? name - '_R1' : name - '_R2'
    """
    demultiplex.py -b ${barcode_file} -f . -s ${params.suffix1} ${params.suffix2} > ${prefix}.log
    """
}