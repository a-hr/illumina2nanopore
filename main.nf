#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include modules
include { basecall } from './modules/basecall'
include { fastqc as raw_fastqc; fastqc as trimmed_fastqc} from './modules/fastqc'
include { merge_fastqs } from './modules/merge_fastqs'
include { multiqc } from './modules/multiqc'
include { stats as adapter_stats; stats as bc_stats; reverse_complement} from './modules/seqkit'
include { adapter_trim; demultiplex_bc } from './modules/demultiplexing'

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
    .value(params.csv_file)
    .set{csv_file}

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

    adapter_demultiplex_logs = adapter_trim.out.log

    // trimmed fastq QC
    adapter_trim.out.trimmed_fastqs \
    | trimmed_fastqc 

    adapter_stats(
        adapter_trim.out.trimmed_fastqs,
        "adapter_trimmed"
    )

    fastqc_logs = fastqc_logs
        .concat(trimmed_fastqc.out)
        .collect()  // makes multiqc wait for all to be emitted

    // reverse complement of trimmed fastqs, filter out files containing unknown in the name
    adapter_trim.out.trimmed_fastqs | flatten \
    | filter { it =~ /^((?!.*(unknown).*).)*$/ } \
    | reverse_complement

    //* notas: podemos quedarnos solo con forward y seguir adelante, o tomar forward y reverse y hacer un reverse complement de ambos y cruzarlos
    //*     - forward_R1 + reverse_R2 = final_R1
    //*     - forward_R2 + reverse_R1 = final_R2
    //* DE MOMENTO LOS UNIMOS

    //! from here on, we only use on pair of reads (forward and reverse)

    // demutiplex by barcode
    demultiplex_bc(
        reverse_complement.out,
        csv_file
    )

    bc_stats(
        demultiplex_bc.out.fastqs,
        "demultiplexed"
    )

    bc_demultiplex_logs = demultiplex_bc.out.log

    // extract UMI

    // map to reference

    // deduplicate

    // multiqc
    multiqc(
        fastqc_logs,
        adapter_demultiplex_logs,
        multiqc_config
    )

}