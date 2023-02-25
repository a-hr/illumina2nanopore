# progreso

* [X] basecalling con guppy
* [X] merge de los reads
* [X] fastqc del merged fastq
* [X] demultiplex de archivos (fw y rv) en funcion de los adaptadores que tienen
* [X] QC de los archivos demultiplexeados
    * [X] fastqc
    * [X] seqkit stats
* [X] reverse complement de los archivos demultiplexeados
* [ ] demultiplex de archivos (fw y rv) en funcion de sus barcodes
* [ ] UMI extraction
* [ ] mapping
* [ ] deduplication
* [ ] analyses

# todo

* abrir y revisar el fastq unknown (no demultiplexados)

# dudas

* tras sacar una muestra del fastq fw (test_fw.fastq):
    * son fw seguro pq tienen cola polyA
    * son fw pq tenian adaptadores fw
    * tienen barcodes reverse!!! --> ¿está mal el archivo de barcodes de Ane?