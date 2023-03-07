process multiqc {
    publishDir "${params.output_dir}/multiqc", mode: "copy"

    input:
        path fastqc_log
        path demultiplex_orientation_log
        path demultiplex_library_log
        path adapter_trim_log
        path demultiplex_barcode_log
        path star_align_log
        path dedup_log
        path featureCounts_log
        path multiqc_custom_config

    output:
        path "multiqc_report.html"

    """
    multiqc --config $multiqc_custom_config .
    """
}