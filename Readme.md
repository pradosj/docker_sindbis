
# Quick Start

`pradosj/sindbis` is a docker container implementing a pipeline to parse FASTQ reads sequenced from SINDBIS infected cells.

The pipeline takes 2 inputs:

  1. a gzipped FASTQ file with extension `<file>.fastq.gz`. This file should contain 106bp-long-reads from SINDBIS infected cells.
  
  2. a tab separated file `<file>.index.tsv` describing how to demultiplex sequences into the dissected areas. The first column must match the pattern `pupName_dissectionName`, and the `dissectionName` of the injection site must contains the string `Inj`. 
  

Assuming the 2 input files exist, the pipeline can be run using the following command line into a terminal:
```
docker run --rm -v $(pwd):/export pradosj/sindbis <file>.sindbis/all
```

## Example `.index.tsv`
Here is an example of input file for `.index.tsv`
```
name	sequence
Q_S1_Inj	ATCACG
Q_S1c	CGATGT
Q_M1	TTAGGC
Q_M1c	TGACCA
Q_S2	ACAGTG
Q_S2c	GCCAAT
Q_Str	CAGATC
Q_Strc ACTTGA
```


# Installation

To run the pipeline, you need to install [Docker](https://www.docker.com/get-started), and make sure enough memory resources are allocated to Docker (mine is set to 12Gb). To check the memory allocation: launch docker and in the docker menu navigate to `Preferences... >> Resources >> ADVANCED`.


# Parameters

You can set the desired filtering parameter using parameter `FILTER=10`. This parameter control the minimum amount of barcode required in the Injections sites (default 100).
```
docker run --rm -v $(pwd):/export pradosj/sindbis FILTER=10 <file>.sindbis/all
```


# Pipeline output

The pipeline is implemented using Makefile rules, and generate multiple files in the ".sindbis/" subdirectory

 1. `demuxmap.tsv`: a table mapping the demultiplexing index sequence to a MOUSE_DISSECTION name. The index is a 6bp nucleotide sequence that will be matched against the reads sequence at position 89-94 for demultiplexing.
 2. `*.fastq.gz`: demultiplexed fastq files
 3. `demux_report.txt`: summary statistics with the number of read in each demultiplexed region
 4. `dissections.txt`: the list of dissection names extracted by parsing `demux_report.txt`
 5. `links.txt`: pairs of Injection-Target dissections to link
 6. `*.umi.class.tsv.gz`: a sorted table of UMI corrected barcodes sequence classified by type. Each row in the table give the number of time the barcode is seen with the same UMI sequence. The `umi` is the 12bp sequence extracted from the reads starting at position 85. `bc32` is a 32bp barcode sequence found in the 32 first bp of the reads. `class` is a classification of the reads according to the content of ` bc32` sequence: if the last 8bp sequence of the barcode match with the sequence `ATCAGTCA` with a tolerance of 1 mismatch, the read is classified as `spike`, otherwise if the last 2bp are `YY` it is classified as `viral`. 
 7. `*.umi.class.viral.fasta`: the set of viral barcodes (bc32) found in `*.umi.class.tsv.gz`. The header of the fasta sequence also contains the umi-corrected count of barcode (i.e. the number of time the barcode is seen with distinct umi sequence).
 8. `MOS3P21_MO_Inj.umi.class.viral.self.bam`: the sequence of above `.fasta` file mapped onto itself with `bowtie` (allowing for 2 mismatches).
 9. `*.bam.clusters.tsv`: result of error correction on the barcode sequences. 
 10. `*.bam.clusters.filter10.fasta`: the error corrected list of barcodes found in the injection site, with at list 10 counts.
 11. `*.self.bam.clusters.filter10.bowtie_aln/`: a directory with the result of mapping the barcodes in the targets dissections onto the barcodes found into the injection site.
     a. `*.self.bam.clusters.filter10.bowtie_aln/*.bam`: result of mapping the barcodes of a given target onto the barcode of the injection site (mapping done with `bowtie`, with a tolerance of 2 mismatches).
     b. `*.self.bam.clusters.filter10.bowtie_aln/*.bam.tsv.gz`: identical to `../*.umi.class.tsv.gz` with the addition of 2 columns where barcodes as been mapped onto barcodes of the injection site.
 


# Generate Barcode Count Matrix for R

For R/Bioconductor users, we have introduce a command to generate a SummarizedExperiment object containing umi-corrected barcode counts. This command looks recursively into the given directory for all files matching pattern `*.umi.class.viral.bam.tsv.gz`; read them to extract the number of barcode into each target that match a barcode of the injection site.
```
docker run --rm -v $(pwd):/export pradosj/sindbis <directory>/umi_corrected_count_matrix.rds
```




