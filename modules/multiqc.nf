process multiqc {
    publishDir "${params.output_dir}/multiqc", mode: "copy"

    input:
        path fastqc_log
        path demultiplex_log
        path multiqc_custom_config
    output:
        path "multiqc_report.html"

    """
    multiqc --config $multiqc_custom_config .
    """
}