#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ------------ MODULES ------------

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
    dedup_UMI;
    cluster_UMI } from './modules/umi_tools'

include { 
    align_STAR;
    align_STAR as consensus_STAR } from './modules/star'

include { 
    align_MINIMAP;
    align_MINIMAP as consensus_MINIMAP } from './modules/minimap'

include { BAM_INDEX } from './modules/samtools'

include { 
    featureCounts as dup_featureCounts;
    featureCounts as dedup_featureCounts } from './modules/subread'

include { plot_results; group_results } from './modules/visualization'


// ------------ FUNCTIONS ------------

def getPrefix(file) {
    return (file.name =~ /(.*)_[0-9]{1,10}\.fast5/)[0][1]
}


// ------------ LOGGING ------------

// copy the files in params.csv_dir to the output directory
process input_logger {
    publishDir "${params.output_dir}/input_logs", mode: 'copy'

    input:
        val csv_dir
        path input_params
        path fastq_dir  // can be either fastq_dir or fast5_dir

    output:
        path "*", includeInputs: true

    script:
        """
        cp ${csv_dir.getParent()}/* .
        ls $fastq_dir > input_files.txt
        """
}


// ------------ INPUT FILES ------------

if (params.enable_basecalling) {
    Channel
        .fromPath("${params.fast5_dir}/*.fast5")
        .set{fast5_files}
    Channel
        .fromPath(params.fast5_dir, type: 'dir')
        .set{fast5_dir}
}

if (!params.enable_isoform_counting) {
    Channel
        .fromPath(params.gtf_file)
        .set{gtf_file}
} else {
    Channel
        .value(1)
        .set{gtf_file}  // dummy file
}


// ------------ RESOURCE FILES ------------

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
    .fromPath(params.lib_csv)
    .set{csv_dir}
Channel
    .value(params.ref_fasta)
    .set{ref_fasta}

// ------------ CONFIG FILES ------------

Channel
    .fromPath(params.multiqc_config)
    .set{multiqc_config}


// ------------ INPUT PARAMETERS ------------

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
    /* ----- LOGGING ----- */
    input_logger(
        csv_dir,
        Channel.fromPath("$baseDir/*.yaml"),
        {params.enable_basecalling ? params.fast5_dir : params.fastq_dir}
    )

    /* ----- BASECALLING ----- */
    if (params.enable_basecalling) {
        prefix = fast5_files.map{ it -> getPrefix(it)}.first()

        fast5_dir \
        | basecall \
        | collect \
        | merge_fastqs \
        | set { fastq_file }
    } else {
        Channel
            .fromPath("$params.fastq_dir/*.fastq.gz") \
            | set { fastq_file }
    }
    
    // raw read QC
    raw_fastqc(fastq_file)
    multiqc_logs = raw_fastqc.out
    
    /* ----- ORIENTATION DEMULTIPLEXING ----- */
    if (params.enable_orientation_demultiplexing) {

        demultiplex_orientation(
            fw_adapters,
            rv_adapters,
            fastq_file
        )

        demultiplex_orientation.out.trimmed_fastqs \
        | oriented_fastqc
        
        orientation_stats(
            demultiplex_orientation.out.trimmed_fastqs,
            "orientation_dmplexed"
        )

        multiqc_logs = multiqc_logs
            .concat(oriented_fastqc.out)
            .collect()  // makes multiqc wait for all to be emitted

        // output channels
        oriented_fastqs     = demultiplex_orientation.out.trimmed_fastqs
        multiqc_logs        = multiqc_logs.concat(demultiplex_orientation.out.log)

    }
    else {
        // output channels
        oriented_fastqs     = fastq_file
    }

    /* ----- REVERSE COMPLEMENT OF TRIMMED READS ----- */
    oriented_fastqs \
    | flatten \
    | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
    | reverse_complement

    reverse_complement.out \
    | flatten \
    | filter { it.name =~ /^(forward_R1|reverse_R2)\.fastq.gz$/ } \
    | collect \
    | merge_forward_reverse \
    | set { complemented_fastqs }
    
    /* LIBRARY DEMULTIPLEXING */

    if (params.enable_library_demultiplexing) {

        demultiplex_library(
            complemented_fastqs,
            lib_csv
        )

        lib_stats(
            demultiplex_library.out.fastqs.collect(),
            "library_dmplexed"
        )

        // output channels
        demultiplex_library.out.fastqs \
        | flatten \
        | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
        | set {lib_dmplexed_fastqs}

        multiqc_logs = multiqc_logs.concat(demultiplex_library.out.log.collect())

    }
    else {
        // output channels
        lib_dmplexed_fastqs = complemented_fastqs
    }


    /* INTERNAL ADAPTER TRIMMING */
    adapter_trim(
        lib_dmplexed_fastqs,
        internal_adapters
    )

    trim_stats(
        adapter_trim.out.fastqs.collect(),
        "adapter_trimmed"
    )

    multiqc_logs = multiqc_logs.concat(adapter_trim.out.log.collect())
    
    /* BARCODE DEMULTIPLEXING */
    if (params.enable_barcode_demultiplexing) {

        demultiplex_bc(
            adapter_trim.out.fastqs,
            bc_csv
        )

        bc_stats(
            demultiplex_bc.out.fastqs.collect(),
            "bc_dmplexed"
        )

        // output channels
        bc_demultiplex_fastqs   = demultiplex_bc.out.fastqs
        multiqc_logs            = multiqc_logs.concat(demultiplex_bc.out.log)

    }
    else {
        // output channels
        bc_demultiplex_fastqs   = adapter_trim.out.fastqs
    }

    /* EXTRACT UMI */

    if (params.enable_UMI_treatment || params.enable_UMI_clustering) {
        bc_demultiplex_fastqs \
        | flatten \
        | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
        | extract_UMI

        // output channels
        extracted_fastqs = extract_UMI.out
    }
    else {
        // output channels
        bc_demultiplex_fastqs \
        | flatten \
        | filter { it.simpleName =~ /^((?!.*(unknown).*).)*$/ } \
        | set { extracted_fastqs }
    }

    /* ALIGNMENT */
    if (params.enable_minimap) {
        align_MINIMAP(ref_fasta, extracted_fastqs)

        // output channels
        alignment_bams = align_MINIMAP.out
    }
    else {
        align_STAR(extracted_fastqs, index_dir)

        multiqc_logs   = multiqc_logs.concat(align_STAR.out.logs.collect())

        // output channels
        alignment_bams = align_STAR.out.bams
    }

    /* DEDUPLICATION OF READS */

    if (params.enable_UMI_treatment) {
        alignment_bams \
        | dedup_UMI

        dedup_UMI.out.dedup_bams \
        | BAM_INDEX

        // output channels
        dedup_bams      = dedup_UMI.out.dedup_bams
        multiqc_logs    = multiqc_logs.concat(dedup_UMI.out.logs.collect())
    } 
    else if (params.enable_UMI_clustering) {
        // cluster reads by UMI
        alignment_bams \
        | cluster_UMI

        // realign clustered reads
        if (params.enable_minimap) {
            consensus_MINIMAP(ref_fasta, cluster_UMI.out)

            // output channels
            dedup_bams = consensus_MINIMAP.out
        }
        else {
            consensus_STAR(
                cluster_UMI.out,
                index_dir
            )

            // output channels
            dedup_bams = consensus_STAR.out.bams
        }

        // index clustered reads
        dedup_bams \
        | BAM_INDEX
    }
    else {
        alignment_bams \
        | BAM_INDEX
    }

    /* EXPRESSION QUANTIFICATION */
    dup_featureCounts(
        alignment_bams.collect(),
        gtf_file,
        "dup"
    )

    multiqc_logs = multiqc_logs.concat(dup_featureCounts.out.logs.collect()).collect()

    if (params.enable_UMI_treatment || params.enable_UMI_clustering) {
        dedup_featureCounts(
            dedup_bams.collect(),
            gtf_file,
            "dedup"
        )

        multiqc_logs = multiqc_logs.concat(dedup_featureCounts.out.logs.collect()).collect()
    }

    /* RESULT CLEANING AND PLOTTING */
    if (params.enable_UMI_treatment || params.enable_UMI_clustering) {
        dedup_featureCounts.out.counts.collect() \
        | concat(dup_featureCounts.out.counts.collect()) \
        | plot_results

        dedup_featureCounts.out.counts.collect() \
        | concat(dup_featureCounts.out.counts.collect()) \
        | group_results
    }
    else {
        dup_featureCounts.out.counts.collect() \
        | plot_results

        dup_featureCounts.out.counts.collect() \
        | group_results
    }

    // multiqc TODO: add input channels in modular way
    multiqc(
        multiqc_logs,
        multiqc_config
    )
}