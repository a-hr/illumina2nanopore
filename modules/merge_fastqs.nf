process merge_fastqs {
    publishDir "${params.output_dir}/fastqs/basecalled",
        mode: 'copy',
        enabled: params.publish_basecalled

    input:
        path fastqs
    output:
        path "out/merged.fastq.gz"

    """
    mkdir out
    cat *.fastq.gz > out/merged.fastq.gz
    """
}

process merge_forward_reverse {
    publishDir "${params.output_dir}/fastqs/complemented",
        mode: 'copy',
        enabled: params.publish_merged_complemented

    input:
        path fastqs
    output:
        path "*.fastq.gz"

    """
    cat *.fastq.gz > final_forward.fastq.gz
    """
}