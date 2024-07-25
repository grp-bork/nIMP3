nIMP3
=====

`IMP` (Integrated Meta-omic Pipeline) is a workflow for metagenomic/metatranscriptomic co-assembly originally developed at the [University of Luxembourg](https://git-r3lab.uni.lu/IMP/IMP)*. `nIMP3` is a nextflow-port of [IMP3](https://git-r3lab.uni.lu/IMP/imp3) developed at EMBL Heidelberg in collaboration with the University of Luxembourg and University of Amsterdam for the [NFDI4Microbiota](https://nfdi4microbiota.de/) framework. nIMP3 is powered by the independent [nevermore](https://github.com/cschu/nevermore) workflow component library.

\* Citation:

```
Narayanasamy, S., Jarosz, Y., Muller, E.E.L. et al. IMP: a pipeline for reproducible reference-independent integrated metagenomic and metatranscriptomic analyses. Genome Biol 17, 260 (2016). https://doi.org/10.1186/s13059-016-1116-8
```

![nevermore_workflow](https://raw.githubusercontent.com/grp-bork/nIMP3/main/docs/nevermore.svg)
![nimp3_workflow](https://raw.githubusercontent.com/grp-bork/nIMP3/main/docs/nimp3_workflow.svg)



### Collaborators

#### Bork Group EMBL Heidelberg
Christian Schudoma, Shahriyar Mahdi Robbani, Daniel Podlesny

#### University of Luxembourg
Oskar Hickl, Patrick May

#### University of Amsterdam
Anna Heintz-Buschart


Dependencies
------------

We recommend running `nIMP3` with Docker/Singularity. By default, it makes use of the biocontainers versions of its dependencies (with the exception of bwa/samtools, s. below)

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




