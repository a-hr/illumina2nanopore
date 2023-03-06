# progreso

## BASECALLERS

* duplex basecalling: para secuenciar pair end con un solo flowcell
    * basecaller Dorado (https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revap_14dec2018/duplex-basecalling)
    * guppy_basecaller_duplex (incluido en guppy)

# todo

* parse and clean input names in csv files:
    * remove double spaces
    * replace spaces and underscores with dashes

# error tracking
* queueSize = 15, submitRateLimit = '1sec', pollInterval = '1m': error_count = 21
* queueSize = 30, submitRateLimit = '1sec', pollInterval = '1m': error_count = 8
* queueSize = 30, submitRateLimit = '1sec', pollInterval = '1m': error_count = 8


# dudas

# Input files

* csvs:
    * `barcodes.csv`: contains the barcodes for demultiplexing
        * structure:
            * `Barcode`: barcode sequence
            * `Sequence`: nucleotide sequence of the barcodes found in the reverse reads (5' end)
            * `Sample`: sample name. DO NOT use underscores `_` in the sample name. Use spaces or dashes `-` instead. Use underscores to separate the sample name from the tissue name.
            * `Reverse`: nucleotide sequence of the barcodes found in the forward reads (3' end)
    * `libraries.csv`: contains the *i7* barcodes for demultiplexing by library
        * structure:
            * `Library`: library name. USE A SINGLE WORD. No spaces, dashes or underscores allowed.
            * `Sequence`: i7 forward sequence
            * `Reverse`: i7 reverse sequence
    * `annotations.saf`: tab separated file containing the target sequence locations to be quantified by featureCounts
        * structure:
            * `Chr`: chromosome
            * `Start`: start position
            * `End`: end position
            * `Strand`: strand
            * `Name`: gene name
            * `Length`: length of the gene

# Options

* Reads per FASTQ file (-q or --records_per_fastq): The number of reads to put in a single FASTQ file (see output format below). Set this to zero to output all reads into one file (per run id, per caller). The default value is 4000.
* 