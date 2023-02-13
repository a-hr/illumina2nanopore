#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process basecall {
    tag "${fast5}"
    publishDir "${params.output_dir}/fastqs", mode: 'copy', saveAs: { out_name }
    
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
    tag "${fast5_dir}"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'
    
    input:
        path fast5_dir
    output:
        path "*_fastqc.{zip,html}"
    
    """
    fastqc --nano -t 8 -o . ${fast5_dir}
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

// Function definitions
def getPrefix(file) {
    return (file.name =~ /(.*)_[0-9]{1,10}\.fast5/)[0][1]
}

// Define variables
Channel.fromPath("${params.fast5_dir}/*.fast5").set{fast5_files}
Channel.fromPath(params.fast5_dir, type: 'dir').set{fast5_dir}

prefix = fast5_files.map{ it -> getPrefix(it)}.first()

workflow {
    
    fast5_files \
    | basecall
    merge_fastqs(prefix, basecall.out.collect())
    
    /*
    fast5_dir \
    | fastqc
    */
}