nevermoreIMP
============

`IMP` (Integrated Meta-omic Pipeline) is a workflow for metagenomic/metatranscriptomic co-assembly originally developed at the [University of Luxembourg](https://git-r3lab.uni.lu/IMP/IMP). `nevermoreIMP` is a nextflow-port of [IMP3](https://git-r3lab.uni.lu/IMP/imp3) powered by the independent [nevermore](https://github.com/cschu/nevermore) workflow component library.


Dependencies
------------

`nevermoreIMP` best runs with Docker/Singularity and, by default, makes use of the biocontainers versions of its dependencies (with the exception of bwa/samtools, s. below)

### Essential/Mandatory

* megahit
* bwa
* samtools

NOTE: `bwa` and `samtools` need to be run from a shared container/environment.

### Optional

* bbmap (bbduk, reformat)
* kraken2
* sortmeRNA
* FastQC
* MultiQC




