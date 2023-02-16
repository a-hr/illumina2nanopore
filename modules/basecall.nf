process basecall {
    tag "${fast5}"
    
    input:
        path fast5
    output:
        path "pass/*.gz"

    script:
        out_name = fast5.simpleName + ".fastq.gz"
        """
        guppy_basecaller -i . -s . -c ${params.config_file} --num_callers 1 --cpu_threads_per_caller 8 --compress_fastq
        """
}