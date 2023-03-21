# progreso

5’_AATGATACGGCGACCACCGAGATCTACACGTTCAG-AGTTCTACAGTCCGACGATC-TARGET-TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC-i7index-ATCTCGTATGCCGTCTTCTGCTTG_3’

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