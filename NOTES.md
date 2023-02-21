# to run

* required config file name: `dna_r9.4.1_450bps_hac` (can get flowcell and kit from html report)
    * To list all supported flowcell and library kits type: ``guppy_basecaller --print_workflows``

* Run command: 

    ```bash
    guppy_basecaller \
        -i inputs/fast5s/ \
        -s fastqs/ \
        -c dna_r9.4.1_450bps_hac.cfg \
        --cpu_threads_per_caller 12 \
        --num_callers 2 \
        --compress_fastq
    ```