#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process basecall {
    tag "${fast5}"
    
    input:
        path fast5
    output:
        path "pass/*.gz"

    script:
        out_name = fast5.simpleName + ".fastq.gz"
        """
        guppy_basecaller -i . -s . -c ${params.config_file} --num_callers 1 --cpu_threads_per_caller 8 --compress_fastq
        """
}

process fastqc {
    tag "${fastq_file}"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    
    input:
        path fastq_file
    output:
        path "*_fastqc.{zip,html}"
    
    """
    fastqc -t 8 -o . ${fastq_file}
    """
}

process merge_fastqs {
    publishDir "${params.output_dir}/fastq", mode: "copy"

    input:
        val sample
        path "${sample}_??.fastq.gz"
    output:
        path "merged_*.fastq.gz"

    """
    cat ${sample}*.fastq.gz > merged_${sample}.fastq.gz
    """
}

process multiqc {
    publishDir "${params.output_dir}/multiqc", mode: "copy"

    input:
        path fastqc_log
        path demultiplex_log
    output:
        path "multiqc_report.html"

    """
    multiqc .
    """
}

process seqkit {
    tag "${fastq_file}"
    publishDir "${params.output_dir}/seqkit", mode: 'copy'
    
    input:
        path fastq_file
    output:
        path "*_seqkit.{txt,html}"
    
    """
    seqkit stats -a -T ${fastq_file} > ${fastq_file.simpleName}_seqkit.txt
    """
}

process adapter_trim {
    // trims the adapters, filters by length and demultiplexes by strand (forward/reverse)
    tag "${fastq_file}"
    publishDir "${params.output_dir}/trimmed_fastqs", pattern: "trimmed_*.fastq.gz", mode: 'copy'
    
    input:
        path fastq_file
    output:
        path "trimmed_*.fastq.gz", emit: trimmed_fastqs
        path "*_demultiplex.log", emit: log
    
    """
    cutadapt \\
        -j 0 \\
        -a forward="AATGATACGGCGACCACCGAGATCTACA...TATGCCGTCTTCTGCTTG" \\
        -a reverse="CAAGCAGAAGACGGCATA...TGTAGATCTCGGTGGTCGCCGTATCATT" \\
        -e 0.2 \\
        -m 50 \\
        -M 300 \\
        -o trimmed_{name}.fastq.gz \\
        ${fastq_file} > ${fastq_file.simpleName}_demultiplex.log
    """
}

// Function definitions
def getPrefix(file) {
    return (file.name =~ /(.*)_[0-9]{1,10}\.fast5/)[0][1]
}

// Define variables
Channel.fromPath("${params.fast5_dir}/*.fast5").set{fast5_files}
Channel.fromPath(params.fast5_dir, type: 'dir').set{fast5_dir}
Channel.fromPath("${params.fastq_dir}/merged_*.fastq.gz").set{fastq_file}

prefix = fast5_files.map{ it -> getPrefix(it)}.first()

workflow {
    
    // fast5_files \
    // | basecall
    // merge_fastqs(prefix, basecall.out.collect())
    
    fastq_file \
    | adapter_trim

    adapter_trim.out.trimmed_fastqs \
    | fastqc

    adapter_trim.out.trimmed_fastqs \
    | seqkit

    fastqc_logs = fastqc.out
    demultiplex_logs = adapter_trim.out.log

    multiqc(fastqc_logs, demultiplex_logs)

}