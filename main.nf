#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include modules
include { basecall } from './modules/basecall'

include { 
    fastqc as raw_fastqc;
    fastqc as trimmed_fastqc} from './modules/fastqc'

include { 
    merge_fastqs;
    merge_forward_reverse} from './modules/merge_fastqs'

include { multiqc } from './modules/multiqc'

include {
    stats as orientation_stats;
    stats as lib_stats;
    stats as bc_stats;
    reverse_complement} from './modules/seqkit'

include { 
    demultiplex_orientation;
    demultiplex_bc;
    demultiplex_library } from './modules/cutadapt'

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
Channel
    .value(params.bc_csv)
    .set{bc_csv}
Channel
    .value(params.lib_csv)
    .set{lib_csv}
Channel
    .value(params.multiqc_config)
    .set{multiqc_config}

// config files
Channel
    .fromPath(params.multiqc_config)
    .set{multiqc_config}

// input parameters
Channel
    .from([params.fw_fp_P5_1, params.fw_tp_P3_1])
    .collect()
    .set{fw_adapters}
Channel
    .from([params.rv_fp_P3_1, params.rv_tp_P5_1])
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
    
    /* 
    - P5 and P3 adapter trimming
    - orientation guess (forward or reverse) based demultiplexing
    */
    demultiplex_orientation(
        fw_adapters,
        rv_adapters,
        fastq_file
    )

    dmplx_orient_logs = demultiplex_orientation.out.log

    /* QC on demultiplexed files */
    demultiplex_orientation.out.trimmed_fastqs \
    | trimmed_fastqc
    
    orientation_stats(
        demultiplex_orientation.out.trimmed_fastqs,
        "orientation_dmplexed"
    )

    fastqc_logs = fastqc_logs
        .concat(trimmed_fastqc.out)
        .collect()  // makes multiqc wait for all to be emitted

    /* REVERSE COMPLEMENT OF TRIMMED READS
    - merge fw_R1 and rv_R2 into final_R1
    - discard fw_R2, rv_R1 and unknown
    */

    demultiplex_orientation.out.trimmed_fastqs \
    | flatten \
    | filter { it =~ /^((?!.*(unknown).*).)*$/ } \
    | reverse_complement

    reverse_complement.out \
    | flatten \
    | filter { it.name =~ /^(forward_R1|reverse_R2)\.fastq.gz$/ } \
    | collect \
    | merge_forward_reverse
    
    /* LIBRARY DEMULTIPLEXING */
    demultiplex_library(
        merge_forward_reverse.out,
        lib_csv
    )

    lib_stats(
        demultiplex_library.out.fastqs,
        "library_dmplexed"
    )

    dmplx_lib_logs = demultiplex_library.out.log

    /* INTERNAL ADAPTER TRIMMING */
    
    /* BARCODE DEMULTIPLEXING */
    // demultiplex_bc(
    //     reverse_complement.out,
    //     bc_csv
    // )

    // bc_stats(
    //     demultiplex_bc.out.fastqs,
    //     "demultiplexed"
    // )

    // bc_demultiplex_logs = demultiplex_bc.out.log

    // extract UMI

    // map to reference

    // deduplicate

    // multiqc
    multiqc(
        fastqc_logs,
        dmplx_orient_logs,
        multiqc_config
    )

}