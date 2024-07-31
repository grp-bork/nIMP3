# nIMP3 workflow
<table>
  <tr width="100%">
    <td width="150px">
      <a href="https://www.bork.embl.de/"><img src="https://www.bork.embl.de/assets/img/normal_version.png" alt="Bork Group Logo" width="150px" height="auto"></a>
    </td>
    <td width="425px" align="center">
      <b>Developed by the <a href="https://www.bork.embl.de/">Bork Group</a> in collaboration with the <a href="https://git-r3lab.uni.lu/IMP/IMP">IMP developers</a></b><br>
      Raise an <a href="https://github.com/grp-bork/nIMP3/issues">issue</a> or <a href="mailto:N4M@embl.de">contact us</a><br><br>
      See our <a href="https://www.bork.embl.de/services.html">other Software & Services</a>
    </td>
    <td width="250px">
      Contributors:<br>
      <ul>
        <li>
          <a href="https://github.com/cschu/">Christian Schudoma</a> <a href="https://orcid.org/0000-0003-1157-1354"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://github.com/mahdi-robbani/">Mahdi Robbani</a> <a href="https://orcid.org/0000-0003-0161-0559"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://github.com/danielpodlesny/">Daniel Podlesny</a> <a href="https://orcid.org/0000-0002-5685-0915"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
      </ul>
    </td>
    <td width="250px">
      Collaborators:<be>
      <ul>
        <li>
          <a href="https://github.com/a-h-b/">Anna Heintz-Buschart</a> <a href="https://orcid.org/0000-0002-9780-1933"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://www.uni.lu/lcsb-en/people/patrick-may/">Patrick May</a> <a href="https://orcid.org/0000-0001-8698-3770"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
        <li>
          <a href="https://www.uni.lu/lcsb-en/people/oskar-hickl/">Oskar Hickl</a> <a href="https://orcid.org/0000-0001-9959-8767"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
        </li>
      </ul>
    </td>
  </tr>
  <tr>
    <td colspan="4" align="center">The development of this workflow was supported by <a href="https://www.nfdi4microbiota.de/">NFDI4Microbiota <img src="https://github.com/user-attachments/assets/1e78f65e-9828-46c0-834c-0ed12ca9d5ed" alt="NFDI4Microbiota icon" width="20px" height="20px"></a> 
</td>
  </tr>
</table>

---
#### Description

`IMP` (Integrated Meta-omic Pipeline) is a workflow for metagenomic/metatranscriptomic co-assembly originally developed at the [University of Luxembourg](https://git-r3lab.uni.lu/IMP/IMP)*. `nIMP3` is a nextflow-port of [IMP3](https://git-r3lab.uni.lu/IMP/imp3) developed at EMBL Heidelberg in collaboration with the University of Luxembourg and University of Amsterdam for the [NFDI4Microbiota](https://nfdi4microbiota.de/) framework. `nIMP3` is powered by the independent [nevermore](https://github.com/cschu/nevermore) workflow component library.

#### Citation
This workflow: badge tbd.

Also cite:
```
Narayanasamy S, Jarosz Y, Muller EE, et al. IMP: a pipeline for reproducible reference-independent integrated metagenomic and metatranscriptomic analyses. Genome Biol. 2016;17(1):260. Published 2016 Dec 16. doi:10.1186/s13059-016-1116-8
```
---
# Overview
![nevermore_workflow](https://raw.githubusercontent.com/grp-bork/nIMP3/main/docs/nevermore.svg)
![nimp3_workflow](https://raw.githubusercontent.com/grp-bork/nIMP3/main/docs/nimp3_workflow.svg)

---
# Requirements
We recommend running `nIMP3` with Docker/Singularity. By default, it makes use of the biocontainers versions of its dependencies (with the exception of `bwa`/`samtools`, see below)

## Essential/Mandatory

* `megahit`
* `bwa`
* `samtools`

NOTE: `bwa` and `samtools` need to be run from a shared container/environment.

## Optional

* `bbmap` (`bbduk`, `reformat`)
* `kraken2`
* `sortmeRNA`
* `FastQC`
* `MultiQC`

---
# Usage
## Cloud-based Workflow Manager (CloWM)
This workflow will be available on the CloWM platform (coming soon).

## Command-Line Interface (CLI)
## Input files
Fastq files are supported and can be either uncompressed (but shouldn't be!) or compressed with `gzip` or `bzip2`. Sample data must be arranged in one directory per sample.

### Per-sample input directories
All files in a sample directory will be associated with the name of the sample folder. Paired-end mate files need to have matching prefixes. Mates 1 and 2 can be specified with suffixes `_[12]`, `_R[12]`, `.[12]`, `.R[12]`. Lane IDs or other read id modifiers have to precede the mate identifier. Files with names not containing either of those patterns will be assigned to be single-ended. Samples consisting of both single and paired end files are assumed to be paired end with all single end files being orphans (quality control survivors). 




