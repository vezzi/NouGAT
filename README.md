# NouGAT
De novo assembly is certainly one of the most difficult tasks of today genomics. On the one hand new sequencing technologies have given the ability to the research community to sequence an increasing number of genome, however, on the other hand, the high number of tools, the difficulties in evaluating the results, and the need of high computational resources makes de novo assembly of genomes still an holy grail.

NouGAT, NGI open universal Genome Assembly Toolbox, is a pipeline that allows, to a certain extent, to automate some of the most complex common processes that take place during the analysis of a de novo assembly project. The pipeline aims to generate a first draft assembly that can be used as a first step towards the production of a better assembly or in order to draw the first biological conclusion.

The pipeline is structured in three different sub-pipelines and in a number of commodity scripts that allow to speed up the analysis in the Uppmax environment.

## Installation
NouGAT is a python package. From a practical point of view it is a wrapper around a certain number of tools that allow to perform de novo assembly analysis.
In order to install it the following steps must be followed:
* Install conda (http://docs.continuum.io/anaconda/install.html#linux-install)
* create a virtual enviorment for the de novo pipeline `conda create -n DeNovoPipeline anaconda`
* activate the virtual enviorment `source activate DeNovoPipeline`
* install report lab package to generate reports `pip install reportlab`
* clone the repository `git clone https://github.com/SciLifeLab/NouGAT.git`
* install de novo pipeline `python setup.py develop`

NouGAT is a wrapper around several tools. The tools needs to be installed and it is highly recommended to have them available on the path.
The supported tools are:
* fastqc
* trimmomatic
* bwa
* picard-tools
* samtools
* ABySS
* ALLPATHS-LG
* cabog
* MaSuRCA
* SOAPdenovo
* SPAdes
* Trinity
* FRC_align --> https://github.com/vezzi/FRC_align
* qaTools   --> https://github.com/vezzi/qaTools

N.B. Not all tools needs to be installed, only the tools that one plans to use. In order to use new tools a new wrapper function needs to be inserted.


## Configuration Files
The pipeline needs two configuration files to be run:
* global configuration: this configuration file contains a description of the sub-pipelines and links to the tools. The repository contains predefined global configurations for milou, nestor, amanita, and picea (n.b. some of the path might still point to my home)
* sample configuration: this describes the samples to be assembled/analysed

### Global Configuration file
This file has two main sections: “Pipelines” and “Tools”.
The former section describes the pipelines implemented and lists which tools can be used in the pipeline.
If pipeline A can use tools T1, T2, and T3 then only these tools can be specified when calling pipeline A, moreover tools T1, T2, and T3 need to be properly installed.

```
 Pipelines:
 QCcontrol: ["fastqc", "abyss", "trimmomatic", "align"]
```

"Tools" section contains an entry for each tool, for each tool the `bin` field must contain the path to the tool (n.b., some tools require the directory where the tool is installed, other require the binary, other require that all the commands are correctly present in the path). `options` is currently not used, however future versions of the pipeline will tak advantage of this field.

```
Tools:
 fastqc:
  bin: /sw/apps/bioinfo/fastqc/0.10.1/milou/fastqc
  options: [--threads,  "16" ,  --outdir,  fastqc]
```

### Sample Configuration file
Sample configuration file specifies which pipeline need to be run (one of those present in the global configuration file), which tools, and in which order. It must be kept in mind that the order in which tools are run can deeply changes the results. Tools cannot be run more than once (i.e., if a tool is present more than once the tool will be run only the first time).

The sample configuration files contains a number of fields (not all mandatory) to be specified, plus a section to describe the library(ies) to be analysed.

More in details:

```
pipeline:
 PIPELINE_TO_BE_RUN ## This field specifies which pipeline needs to be run. Only implemented pipelines can be run.
tools:
 [T1,T2,T3,...] ## Which tools need to be used, and in which order.
output: OUTPUT_NAME  #prefix to append to all output files
minCtgLength: MIN_CTG_LGTH ## minimum contig length to be considered. This parameter is used in several steps in order to discard short contigs (default 2000)
genomeSize: EXP_GENOME_SIZE ## expected genome size (rough estimation)
threads : NUM_THREADS ## number of threads to be used in parallel steps
kmer: KMER ## in case a too needs a predefined kmer used this (suggested is 54, but it is a rule-of-thumb)
reference: PATH_TO_REFERENCE ## this field is mandatory when an alignment needs to be executed, depending on the pipeline it can be the de novo assembly to be evaluated or a reference genome (in case of Mate Pairs)
```


The next section of sample configuration file is 'libraries' and contains a description of the libraries and of the sequencing runs:
libraries:

Each entry (lib1, lib2, etc.) contains the following mandatory fields:
 
```
libN:
 pair1: PATH_TO_PAIR_1 ## Path to first pair
 pair2: PATH_TO_PAIR_2 ## Path to second pair
 orientation: PAIR_ORIENTATION ## Pair read orientation (innie or outtie)
 insert: INSERT_SIZE ## insert size (expected)
 std: STANDARD_DEVIATION ## standard deviation of the insert size  (expected)
```
There is also support for single ended libraries.

It is important to note that (despite the name) lib1, lib2, … identify different sequencing runs. In reality the concept of “library” is represented by the insert size, i.e., lib entries libi and libj with i not equal to j are considered to be part of the same library if and only if the insert size is the same.

### Example
Assemble step, n.b., no reference is specified as it is not needed in this pipeline. This is a typical NGI-S example, of a de novo project (J.Dohe) that is splitted into two projects (J.Dohe_14_01  and J.Dohe_14_02) the former being the Paired end library, and the latter being the Mate pair library. This is often call “allpaths recipe” as this is the assembling strategy suggested by BROAD to use allpaths-lg assembler.

```
pipeline:
 assemble
tools:
 [allpaths]  ## several tools can be specified here
output: J.Dohe
minCtgLength: 2000
genomeSize: 450000000
threads : 16
kmer: 54
reference: 
libraries:
 lib1:
  pair1: /proj/a2010002/INBOX/J.Dohe_14_01/P101/HISEQ_RUN/read_1.fastq.gz
  pair2: /proj/a2010002/INBOX/J.Dohe_14_01/P101/HISEQ_RUN/read_2.fastq.gz
  orientation: innie
  insert: 180
  std: 30
 lib2:
  pair1: /proj/a2010002/INBOX/J.Dohe_14_02/P101/HISEQ_RUN/read_1.fastq.gz
  pair2: /proj/a2010002/INBOX/J.Dohe_14_02/P101/HISEQ_RUN/read_2.fastq.gz
  orientation: outtie
  insert: 3000
  std: 500
```


## The Pipelines
De_Novo_Scilife currently implements 3 different pipelines:

* QCcontrol
* assemble
* evaluate

there is also a fourth pipeline (align) that is a sort of dummy pipeline to align reads against a reference. The align pipeline can be run as a part of other pipelines and makes easyer the specification of pipelines.

Pipelines are described in the global configuration file, in particular this file tells which tools can be run when employing a pipeline. The tools that will be effectively used and the order is specified by the sample configuration file.

Currently the pipelines implement the following tools:
```
 QCcontrol: ["fastqc", "abyss", "trimmomatic", "align"]
 assemble: ["abyss", "allpaths", "cabog",  "masurca", "soapdenovo", "spades", "trinity"]
 evaluete: ["align", "qaTools", "FRC"]
 align: []
```

### QCcontrol
This pipeline implements the instruments needed to perform Quality Control with de novo samples. In particular, it allows access to 4 instruments

``` QCcontrol: ["fastqc", "abyss", "trimmomatic", "align"] ```


A typical way to specify this pipeine when processing a sample is the following:

```
 pipeline:
  QCcontrol
 tools:
   ['trimmomatic', 'fastqc', 'abyss']
```

this will trim the reads, perform fastqc on the trimmed reads and plot a kmer graph using abyss utility for it.
In case a reference is available also “align” can be specified. 

### assemble

### evaluate

### How to run: example





