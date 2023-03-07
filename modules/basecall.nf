process basecall {
    tag "${fast5}"

    input:
        path fast5_dir
    output:
        path "pass/*.gz"

    """
    guppy_basecaller \\
    -i ${fast5_dir} \\
    -s . \\
    -c ${params.config_file} \\
    --device auto \\
    --compress_fastq \\
    --num_callers ${task.cpus} \\
    --chunks_per_caller 400
    """
}