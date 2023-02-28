# progreso

* [X] basecalling con guppy
* [X] merge de los reads
* [X] fastqc del merged fastq
* [X] demultiplex de archivos (fw y rv) en funcion de los adaptadores que tienen
* [X] QC de los archivos demultiplexeados
    * [X] fastqc
    * [X] seqkit stats
* [X] reverse complement de los archivos demultiplexeados
    * [X] unir fw_R1 y rv_R2 (forward reads)
    * [X] descartar los archivos que no son forward reads (fw_R2 y rv_R1)
* [X] demultiplex de los forward reads en funcion de su libreria
* [ ] trimmado de los adaptadores P5 y P3
* [ ] demultiplex en funcion de barcodes
* [ ] UMI extraction
* [ ] mapping
* [ ] deduplication
* [ ] analyses

# todo



# dudas

* tras sacar una muestra del fastq fw (test_fw.fastq):
    * son fw seguro pq tienen cola polyA
    * son fw pq tenian adaptadores fw
    * tienen barcodes reverse!!! --> ¿está mal el archivo de barcodes de Ane?

# Options