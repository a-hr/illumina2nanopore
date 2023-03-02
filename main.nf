#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include modules
include { basecall } from './modules/basecall'

include { 
    fastqc as raw_fastqc;
    fastqc as oriented_fastqc} from './modules/fastqc'

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
    demultiplex_library;
    adapter_trim } from './modules/cutadapt'

include {
    extract_UMI;
    dedup_UMI } from './modules/umi_tools'

include { STAR_ALIGN } from './modules/star'

include { BAM_INDEX } from './modules/samtools'

include { featureCounts } from './modules/subread'

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
    .value(params.index_dir)
    .set{index_dir}
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
Channel
    .value(params.saf_file)
    .set{saf_file}

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
Channel
    .from([params.fw_P5_2, params.fw_P3_2])
    .collect()
    .set{internal_adapters}

prefix = fast5_files.map{ it -> getPrefix(it)}.first()

workflow {
    /* ----- BASECALLING ----- */
    // fast5_files \
    // | basecall
    // merge_fastqs(prefix, basecall.out.collect())
    
    // raw read QC
    raw_fastqc(fastq_file)
    fastqc_logs = raw_fastqc.out
    
    /* ----- ORIENTATION DEMULTIPLEXING ----- */
    /*
    - demultiplex fw and rv reads based on orientation
    - remove P5 and P3 adapters
    - QC on demultiplexed files
    */
    demultiplex_orientation(
        fw_adapters,
        rv_adapters,
        fastq_file
    )

    dmplx_orient_logs = demultiplex_orientation.out.log

    demultiplex_orientation.out.trimmed_fastqs \
    | oriented_fastqc
    
    orientation_stats(
        demultiplex_orientation.out.trimmed_fastqs,
        "orientation_dmplexed"
    )

    fastqc_logs = fastqc_logs
        .concat(oriented_fastqc.out)
        .collect()  // makes multiqc wait for all to be emitted

    /* ----- REVERSE COMPLEMENT OF TRIMMED READS ----- */
    /*
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

    dmplx_lib_logs = demultiplex_library.out.log.collect()

    // prepare the channel for the next step
    demultiplex_library.out.fastqs \
    | flatten \
    | filter { it.simpleName =~ /^(?!pool_unknown).*$/ } \
    | set {lib_dmplexed_fastqs}

    /* INTERNAL ADAPTER TRIMMING */
    adapter_trim(
        lib_dmplexed_fastqs,
        internal_adapters
    )

    adapter_trim_logs = adapter_trim.out.log.collect()
    
    /* BARCODE DEMULTIPLEXING */
    demultiplex_bc(
        adapter_trim.out.fastqs,
        bc_csv
    )

    bc_demultiplex_logs = demultiplex_bc.out.log

    bc_stats(
        demultiplex_bc.out.fastqs.collect(),
        "bc_dmplexed"
    )

    // extract UMI
    demultiplex_bc.out.fastqs \
    | flatten \
    | extract_UMI

    // read alignment
    STAR_ALIGN(extract_UMI.out, index_dir)

    alignment_multiqc = STAR_ALIGN.out.logs.collect()

    // deduplication of aligned bams
    STAR_ALIGN.out.bams \
    | dedup_UMI
    
    dedup_multiqc = dedup_UMI.out.logs.collect()
    
    // reindexing of bams for IGV
    dedup_UMI.out.dedup_bams \
    | BAM_INDEX

    // expression quantification
    featureCounts(dedup_UMI.out.dedup_bams.collect(), saf_file) 
    featureCounts_multiqc = featureCounts.out.logs.collect()

    // multiqc
    multiqc(
        fastqc_logs,
        dmplx_orient_logs,
        multiqc_config
    )

}