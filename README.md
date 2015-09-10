# NouGAT
De novo assembly is certainly one of the most difficult tasks of today genomics. On the one hand new sequencing technologies have given the ability to the research community to sequence an increasing number of genome, however, on the other hand, the high number of tools, the difficulties in evaluating the results, and the need of high computational resources makes de novo assembly of genomes still an holy grail.

NouGAT, NGI open universal Genome Assembly Toolbox, is a pipeline that allows, to a certain extent, to automate some of the most complex common processes that take place during the analysis of a de novo assembly project. The pipeline aims to generate a first draft assembly that can be used as a first step towards the production of a better assembly or in order to draw the first biological conclusion.

The pipeline is structured in three different sub-pipelines, and a number of convenience [scripts](/sciLifeLab_utils) aimed only for NGI users running their analysis in the UPPMAX environment.


## Table Of Contents
* [Installation](README.md#installation)
 * [Third party tools](README.md#third-party-tools) 
* [Configuration files](README.md#configuration-files)
 * [Global configuration file](README.md#global-configuration-file)
 * [Sample configuration file](README.md#sample-configuration-file)
   * [Example file](README.md#example-file)
* [The pipelines](README.md#the-pipelines)
 * [Quality control](README.md#qccontrol)
 * [Assemble](README.md#assemble)
 * [Evaluate](README.md#evaluate)
 * [NGI specific scripts](README.md#ngi-specific-scripts)
* [How to get started](README.md#how-to-get-started)


## Installation
NouGAT is a python package. From a practical point of view it is a wrapper around a certain number of tools that allow to perform de novo assembly analysis. In order to install NouGAT the following steps must be followed:

```bash
# First we create a new virtual python environment using conda, get it here: 
# http://docs.continuum.io/anaconda/install.html#linux-install
conda create -n DeNovoPipeline anaconda

# Install NouGAT in the new environment
source activate DeNovoPipeline
git clone https://github.com/SciLifeLab/NouGAT.git
cd NouGAT && python setup.py develop
```
**NGI / UPPMAX users only:** See the separate [README](sciLifeLab_utils/README.md) on how to install the NGI specific scripts.

### Third party tools
NouGAT is a wrapper around several tools. The tools needs to be installed and it is highly recommended to have them available on the path. The following table shows the supported tools and their versions:

Tool              | Version                             | Function
----------------- | ----------------------------------- | ----------
FastQC            | 0.11.2                              | Read pre-processing
Trimmomatic       | 0.30                                | Read pre-processing
bwa               | > 0.7.4                             | Read alignment
picard-tools      | < 1.124                             | Read alignment
SAMtools          | 1.1.19                              | Read alignment
KmerGenie         | 1.6741                              | Read pre-processing
ABySS             | 1.3.5                               | Assembler
ALLPATHS-LG       | >= 47918                            | Assembler
CABOG             | 8.1                                 | Assembler
MaSuRCA           | 2.3.2                               | Assembler
SOAPdenovo        | 2.04-r240                           | Assembler
SPAdes            | >= 3.0.0                            | Assembler
Trinity           | >= 2.0.2                            | Assembler
FRC_align         | https://github.com/vezzi/FRC_align  | Assembly evaluation
qaTools           | https://github.com/vezzi/qaTools    | Assembly evaluation

**Note!** Not all tools needs to be installed, only the tools that one plans to use. In order to use tools other than the ones in the list, a new wrapper function needs to be be written.

**Also note!** The versions given are the ones tested by the authors, however others might work if the commandline and/or outputs have not changed significantly (e.g. trinity pre-2.0.2)


## Configuration files
The pipeline needs two configuration files to be run:
* [Global configuration](README.md#global-configuration-file): this configuration file contains a description of the sub-pipelines and links to the tools. The repository contains predefined global configurations for milou, nestor, amanita, and picea (n.b. some of the path might still point to my home)
* [Sample configuration](README.md#sample-configuration-file): this describes the samples to be assembled/analysed

### Global configuration file
This file has two main sections: “Pipelines” and “Tools”. The former section describes the pipelines implemented and lists which tools can be used in the pipeline. If pipeline A can use tools T1, T2, and T3 then only these tools can be specified when calling pipeline A, moreover tools T1, T2, and T3 need to be properly installed. 

"Tools" section contains an entry for each tool, for each tool the `bin` field must contain the path to the tool (n.b., some tools require the directory where the tool is installed, other require the binary, other require that all the commands are correctly present in the path). In `options`, commandline parameters can be given (where supported), see the the documentation for each individual tool.

```yaml
Pipelines:
 QCcontrol: ["fastqc", "abyss", "trimmomatic", "align"]
Tools:
 fastqc:
  bin: /sw/apps/bioinfo/fastqc/0.10.1/milou/fastqc
  options: [--threads,  "16" ,  --outdir,  fastqc]
 abyss:
  ...
```
For a complete example of this file, see the [NGI-specific configuration](config_files/config_de_novo_uppmax.yaml).

### Sample configuration file
Sample configuration file specifies which pipeline need to be run (one of who is present in the global configuration file), which tools are run and in which order. It must be kept in mind that the order in which tools are run can deeply change the results. Tools cannot be run more than once (i.e., if a tool is present more than once the tool will be run only the first time).

The sample configuration files contains a number of fields (not all mandatory) to be specified, plus a section to describe the library(ies) to be analysed. More in details:

```yaml
pipeline:
 # This field specifies which pipeline needs to be run. Only implemented pipelines can be run.
 PIPELINE_TO_BE_RUN
# Which tools need to be used, and in which order.
tools:
 [T1,T2,T3,...]
# Prefix to append to all output files
output: OUTPUT_NAME
# Minimum contig length to be considered. This parameter is used in several steps in order to 
# discard short contigs (default 1000)
minCtgLength: MIN_CTG_LGTH
# Expected genome size. Used by some assemblers, if unsure, give a rough estimate.
genomeSize: EXP_GENOME_SIZE
# Number of threads to be used in parallel steps.
threads : NUM_THREADS
# In case a tool needs a predefined kmer size (Rule-of-thumb, try 61 first)
kmer: KMER
# This field is mandatory when an alignment needs to be executed, e.g, the tools in the evaluation  
# pipeline use the .bam alignments to work. 
reference: PATH_TO_REFERENCE
```

The next section of sample configuration file is 'libraries' and contains a description of the libraries and paths to the fastq files. Each entry (lib1, lib2, etc.) contains the following mandatory fields:
 
```yaml
libN:
 pair1: PATH_TO_PAIR_1 # Path to first pair
 pair2: PATH_TO_PAIR_2 # Path to second pair. Leave this blank in case of single-ended lib.
 orientation: PAIR_ORIENTATION # Pair read orientation (innie or outtie)
 insert: INSERT_SIZE # insert size (expected)
 std: STANDARD_DEVIATION # standard deviation of the insert size  (expected)
```

It is important to note that (despite the name) lib1, lib2, … identify different sequencing runs. In reality the concept of “library” is represented by the insert size, i.e., library entries lib*i* and lib*j* with *i* not equal to *j* are considered to be part of the same library if and only if the insert size is the same.

##### Example file
Sample configuration for the assemble pipeline, note that no reference is specified as it is not needed. This is an example of a typical NGI Stockholm de novo project (J.Dohe) that is split into two projects (J.Dohe_14_01  and J.Dohe_14_02) the former being the paired end, and the latter being the mate pair library. This is often called the “allpaths recipe” as this is the assembly strategy suggested by BROAD for their ALLPATHS-LG assembler.

```yaml
pipeline:
 assemble
tools:
 [allpaths]
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


## The pipelines
NouGAT currently implements 3 different pipelines named:
* `QCcontrol`
* `assemble`
* `evaluate`

There is also a fourth pipeline, `align`, which is a sort of dummy pipeline to align reads against a reference. The align pipeline can be run as a part of other pipeline, this makes it easier to specify custom align pipelines, e.g. one using something other than BWA.

### Quality control
The pipeline is for pre-processing the sequenced reads, i.e. to ascertain the quality of the sequencing experiment and inform decisions to be made for the assembly step. An example would be to run the following tools in succession:

* Trimmomatic, to remove adapter read-through and low quality ends.
* FastQC, to give easy to read quality checks for Q-values, read lengths, etc.
* Abyss. This is for generating k-mer counts from the reads. The pipeline will automatically plot these.

Also useful is to run `align` if a reference assembly for the organism exists. The pipeline will pull useful statistics from the .bam aligned files using samtools and picard, e.g. estimated insert size and duplication rate metrics.

### Assemble
The assembly pipeline will sequentially execute the assembly tools specified in the sample config file. From our earlier [J.Dohe example](README.md#example-file), we in addition want to run ABySS and SOAPdenovo before ALLPATHS-LG:

```yaml
pipeline:
 assemble
tools:
 [abyss, soapdenovo, allpaths]
output: J.Dohe
...
```

### Evaluate
This pipeline is for evaluating the relative quality of the assemblies produced. This includes the standard contiguity metrics (N50, N80, max contig length, etc.) and estimates of the level of mis-assembly. It is recommended to run these tools in the following order:

* Align
* qaTools, to visualise ctg. length vs. coverage vs. GC content, etc.
* FRC, feature response curves to estimate the relative levels of mis-assembly.
 
Note! FRC only accepts two libraries so we recommend to choose one (overlapping) paired-end and one long insert (mate pair) library.

### NGI specific scripts
Described [here](/sciLifeLab_utils).

## How to get started
Coming soon.
