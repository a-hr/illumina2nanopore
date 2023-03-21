process demultiplex_orientation {
    tag "${fastq_file}"

    publishDir "${params.output_dir}/fastqs/oriented",
        pattern: "*.fastq.gz",
        mode: 'copy',
        enabled: params.publish_orient_demultiplexed
    
    input:
        tuple val(five_prime_fw), val(three_prime_fw)
        tuple val(five_prime_rv), val(three_prime_rv)
        path fastq_file
    output:
        path "*.fastq.gz", emit: trimmed_fastqs
        path "*_dmplex.log", emit: log

    """
    cutadapt \\
        -j 0 \\
        -a forward="${five_prime_fw}...${three_prime_fw}" \\
        -a reverse="${five_prime_rv}...${three_prime_rv}" \\
        -e 0.2 \\
        -m ${params.min_length} \\
        -M ${params.max_length} \\
        -o {name}.fastq.gz \\
        ${fastq_file} > ${fastq_file.simpleName}_orient_dmplex.log
    """
}

process demultiplex_library {
    tag { fastq.name }
    publishDir "${params.output_dir}/fastqs/lib_demultiplexed",
        pattern: "*fastq.gz",
        mode: 'copy',
        enabled: params.publish_lib_demultiplexed

    input:
        path fastq
        path lib_csv

    output:
        path "*.gz", emit: fastqs
        path "*.log", emit: log

    script:
    """
    demultiplex_lib.py -b ${lib_csv} -f ${fastq} > ${fastq.simpleName}_lib_dmplex.log
    """
}

process demultiplex_bc {
    tag { fastq.name }

    publishDir "${params.output_dir}/fastqs/bc_demultiplexed",
        mode: 'copy',
        enabled: params.publish_bc_demultiplexed

    input:
        path fastq
        path barcode_file

    output:
        path "*.gz", emit: fastqs
        path "*.log", emit: log

    script:
    prefix = fastq.simpleName
    """
    demultiplex_SE.py -b ${barcode_file} -f . > dmplex_bc_${prefix}.log
    """
}

process adapter_trim {
    tag { fastq.name }
    publishDir "${params.output_dir}/fastqs/trimmed",
        pattern: "*fastq.gz",
        mode: 'copy',
        enabled: params.publish_trimmed

    input:
        path fastq
        tuple val(five_prime), val(three_prime)

    output:
        path "trimmed_*fastq.gz", emit: fastqs
        path "*.log", emit: log

    """
    cutadapt \\
        -j 0 \\
        -e 0.2 \\
        --no-indels \\
        -g ${five_prime}...${three_prime} \\
        -o trimmed_${fastq} \\
        ${fastq} > adapter_trim_${fastq.name}.log
    """
}