
# Quick Start

`pradosj/sindbis` is a docker container implementing a pipeline to parse FASTQ reads sequenced from SINDBIS infected cells.

The pipeline takes 2 inputs:

  1. a gzipped FASTQ file with extension `<file>.fastq.gz`. This file should contain 106bp-long-reads from SINDBIS infected cells.
  
  2. a tab separated file `<file>.index.tsv` describing how to demultiplex sequences into the dissected areas. The first column must match the pattern `pupName_dissectionName`, and the `dissectionName` of the injection site must be equal to `S1InjSite`. 
  

Assuming the 2 input files exist, the pipeline can be run using the following command line into a terminal:
```
docker run --rm -v $(pwd):/export pradosj/sindbis <file>.sindbis/all
```

# Example `.index.tsv`
Here is an example of input file for `.index.tsv`
```
name	primer	sequence
Q_S1InjSite	RT1	ATCACG
Q_S1c	RT2	CGATGT
Q_M1	RT3	 TTAGGC
Q_M1c	RT4	TGACCA
Q_S2	RT5	ACAGTG
Q_S2c	RT6	GCCAAT
Q_Str	RT7	CAGATC
Q_Strc	RT8	ACTTGA
```


# Installation

To run the pipeline, you just need to install [Docker](https://www.docker.com/get-started), and make sure enough memory resources are allocated to Docker (mine is set to 12Gb). To check the memory allocation: launch docker and in the docker menu navigate to `Preferences... >> Resources >> ADVANCED`.

