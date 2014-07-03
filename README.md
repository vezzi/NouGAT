# de_novo_scilife
De novo assembly is certainly one of the most difficult tasks of today genomics. On the one hand new sequencing technologies have given the ability to the research community to sequence an increasing number of genome, however, on the other hand, the high number of tools, the difficulties in evaluating the results, and the need of high computational resources makes de novo assembly of genomes still an holygrail.

De Novo Scilife is a pipeline that allows, to a certain extent, to automate some of the most complex common processes that take place during the analysis of a de novo assembly project. The pipeline aims to generate a first draft assembly that can be used as a first step towards the production of a better assembly or in order to draw the first biological conclusion.

The pipeline is structured in three different sub-pipelines and in a number of commodity scripts that allow to speed up the analysis in the Uppmax environment.

## Installation
De_novo_scilife is a python package. From a practical point of view it is a wrapper around a certain number of tools that allow to perform de novo assembly analysis.
In order to install it the following steps must be followed:
* Install conda (http://docs.continuum.io/anaconda/install.html#linux-install)
* create a virtual enviorment for the de novo pipeline `conda create -n DeNovoPipeline anaconda`
* activate the virtual enviorment `source activate DeNovoPipeline`
* install report lab package to generate reports `pip install reportlab`
* clone the repository `git clone https://github.com/SciLifeLab/de_novo_scilife.git`
* install de novo pipeline `python setup.py develop`

De_novo_pipeline is a wrapper around several tools. The tools needs to be installed and it is highly recommended to have them available on the path.
The supported tools are:
* fastqc
* trimmomatic
* bwa
* picard-tools
* samtools
* abyss
* allpathslg
* cabog
* MaSuRCA
* SOAPdenovo
* spades
* FRC_align --> https://github.com/vezzi/FRC_align
* qaTools   --> https://github.com/vezzi/qaTools

N.B. Not all tools needs to be installed, only the tools that one plans to use. In order to use new tools a new wrapper function needs to be inserted.


## Configuration Files
The pipeline needs two configuration files to be run:
* global configuration: this configuration file contains a description of the sub-pipelines and links to the tools. The repository contains predefined global configurations for milou, nestor, amanita, and picea (n.b. some of the path might still point to my home)
* sample configuration: this describes the samples to be assembled/analysed


### Global Configuration file



### Sample Configuration file
