#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include modules
include { basecall } from './modules/basecall'
include { fastqc as raw_fastqc; fastqc as trimmed_fastqc} from './modules/fastqc'
include { merge_fastqs } from './modules/merge_fastqs'
include { multiqc } from './modules/multiqc'
include { demultiplex_stats; reverse_complement } from './modules/seqkit'
include { adapter_trim } from './modules/trimming'

// function definitions
def getPrefix(file) {
    return (file.name =~ /(.*)_[0-9]{1,10}\.fast5/)[0][1]
}

// ------------ VARIABLES ------------
// paths
Channel
    .fromPath("${params.fast5_dir}/*.fast5")
    .set{fast5_files}
Channel
    .fromPath(params.fast5_dir, type: 'dir')
    .set{fast5_dir}
Channel
    .fromPath("${params.fastq_dir}/merged_*.fastq.gz")
    .set{fastq_file}

// config files
Channel
    .fromPath(params.multiqc_config)
    .set{multiqc_config}

// input parameters
Channel
    .from([params.five_prime_fw, params.three_prime_fw])
    .collect()
    .set{fw_adapters}
Channel
    .from([params.five_prime_rv, params.three_prime_rv])
    .collect()
    .set{rv_adapters}

prefix = fast5_files.map{ it -> getPrefix(it)}.first()

workflow {
    // fast5_files \
    // | basecall
    // merge_fastqs(prefix, basecall.out.collect())
    
    // raw read QC
    raw_fastqc(fastq_file)
    fastqc_logs = raw_fastqc.out
    
    // adapter trimming/demultiplexing
    adapter_trim(
        fw_adapters,
        rv_adapters,
        fastq_file
    )

    demultiplex_logs = adapter_trim.out.log

    // trimmed fastq QC
    adapter_trim.out.trimmed_fastqs \
    | (trimmed_fastqc & demultiplex_stats)

    fastqc_logs = fastqc_logs.concat(trimmed_fastqc.out).collect()  // makes multiqc wait for all to be emitted

    // reverse complement of trimmed fastqs
    adapter_trim.out.trimmed_fastqs | flatten \
    | reverse_complement

    // multiqc
    multiqc(
        fastqc_logs,
        demultiplex_logs,
        multiqc_config
    )

}