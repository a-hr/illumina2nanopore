#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process basecall {
    tag "${fast5}"
    publishDir "${params.output_dir}/fastqs", mode: 'move', saveAs: { out_name }
    
    input:
        path fast5
    output:
        path "pass/*.gz"
    script:
        out_name = fast5.simpleName + ".fastq.gz"
        """
        guppy_basecaller -i . -s . -c dna_r9.4.1_450bps_hac.cfg --num_callers 1 --cpu_threads_per_caller 8 --compress_fastq
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

// Define variables
Channel.fromPath("${params.fast5_dir}/*.fast5").set{fast5_files}
Channel.fromPath(params.fast5_dir, type: 'dir').set{fast5_dir}

ignore_basecall = true

workflow {
    if(!ignore_basecall) {
        fast5_files \
        | basecall
    }

    fast5_dir \
    | fastqc
}