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
    stats as trim_stats;
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

include { 
    featureCounts as fw_featureCounts;
    featureCounts as rv_featureCounts } from './modules/subread'

// function definitions
def getPrefix(file) {
    return (file.name =~ /(.*)_[0-9]{1,10}\.fast5/)[0][1]
}

// ------------ VARIABLES ------------

// paths
if (params.basecall) {
    Channel
        .fromPath("${params.fast5_dir}/*.fast5")
        .set{fast5_files}
    Channel
        .fromPath(params.fast5_dir, type: 'dir')
        .set{fast5_dir}
}
else {
    Channel
        .fromPath("$params.fastq_dir/*.fastq.gz")
        .set {fastq_file}
}

Channel
    .value(params.index_dir)
    .set{index_dir}
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

workflow {
    /* ----- BASECALLING ----- */

    if (params.basecall) {
        prefix = fast5_files.map{ it -> getPrefix(it)}.first()

        fast5_dir \
        | basecall \
        | collect \
        | merge_fastqs \
        | set { fastq_file }
    }
    
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
    | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
    | reverse_complement

    if (params.merge_all_fw) {
        reverse_complement.out \
            | flatten \
            | filter { it.name =~ /^(forward_R1|reverse_R2)\.fastq.gz$/ } \
            | collect \
            | merge_forward_reverse \
            | set { complemented_fastqs }
    }
    else {
        reverse_complement.out \
            | flatten \
            | filter { it.name =~ /^(forward_R1|reverse_R2)\.fastq.gz$/ } \
            | set { complemented_fastqs }
    }
    
    /* LIBRARY DEMULTIPLEXING */
    demultiplex_library(
        complemented_fastqs,
        lib_csv
    )

    lib_stats(
        demultiplex_library.out.fastqs.collect(),
        "library_dmplexed"
    )

    dmplx_lib_logs = demultiplex_library.out.log.collect()

    // prepare the channel for the next step
    demultiplex_library.out.fastqs \
    | flatten \
    | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
    | set {lib_dmplexed_fastqs}

    /* INTERNAL ADAPTER TRIMMING */
    adapter_trim(
        lib_dmplexed_fastqs,
        internal_adapters
    )

    trim_stats(
        adapter_trim.out.fastqs.collect(),
        "adapter_trimmed"
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

    /* EXTRACT UMI */

    demultiplex_bc.out.fastqs \
    | flatten \
    | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
    | extract_UMI

    /* ALIGNMENT */
    STAR_ALIGN(extract_UMI.out, index_dir)

    alignment_multiqc = STAR_ALIGN.out.logs.collect()

    /* DEDUPLICATION OF READS */
    STAR_ALIGN.out.bams \
    | dedup_UMI
    
    dedup_multiqc = dedup_UMI.out.logs.collect()
    
    /* BAM INDEXING */
    dedup_UMI.out.dedup_bams \
    | BAM_INDEX

    /* EXPRESSION QUANTIFICATION */

    // separate bams into forward and reverse
    dedup_UMI.out.dedup_bams \
    | branch {
        forward: it.name =~ /.*forward.*/
        reverse: it.name =~ /.*reverse.*/
    } \
    | set {dedup_filt_bams}

    STAR_ALIGN.out.bams \
    | branch {
        forward: it.name =~ /.*forward.*/
        reverse: it.name =~ /.*reverse.*/
    } \
    | set {star_filt_bams}

    fw_featureCounts(
        dedup_filt_bams.forward.collect(),
        star_filt_bams.forward.collect(),
        saf_file
    )

    rv_featureCounts(
        dedup_filt_bams.reverse.collect(),
        star_filt_bams.reverse.collect(),
        saf_file
    )
    
    featureCounts_multiqc = fw_featureCounts.out.logs.collect()

    // multiqc
    multiqc(
        fastqc_logs,
        dmplx_orient_logs,
        dmplx_lib_logs,
        adapter_trim_logs,
        bc_demultiplex_logs,
        alignment_multiqc,
        dedup_multiqc,
        featureCounts_multiqc,
        multiqc_config
    )

}