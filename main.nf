#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process basecall {
    tag "${fast5}"
    publishDir "${params.output_dir}", mode: 'move'
    
    input:
        path fast5 
    output:
        path "*.gz"
  
    """
    guppy_basecaller -i . -s . -c dna_r9.4.1_450bps_hac.cfg --num_callers 1 --cpu_threads_per_caller 8 --compress_fastq
    """
}

workflow {
    // get the fast5 files from the input directory
    Channel.fromPath("${params.fast5_dir}/*.fast5") \
    | basecall
}