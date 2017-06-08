#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                             NouGat
========================================================================================
 #### Homepage / Documentation
 https://github.com/SciLifeLab/ 
 #### Authors
 Francesco Vezzi
 Remi ...
----------------------------------------------------------------------------------------
*/

revision = grabRevision()

// Pipeline version
version = 0.1

if (!isAllowedParams(params)) {exit 1, "params is unknown, see --help for more information"}


if (params.help) {
  helpMessage(version, revision)
  exit 1
}

if (params.version) {
  versionMessage(version, revision)
  exit 1
}

// Configurable variables
params.align  = false
bwa_index = false
params.genome = false
params.bwa_index =  params.genome ? params.genomes[ params.genome ].bwa ?: false : false
bwa_index_sa = false
bwa_index_bwt = false
bwa_index_pac = false
bwa_index_ann = false
bwa_index_amb = false


if( params.align ) {
    if( !params.bwa_index ){
        exit 1, "No reference genome specified!"
    }
    if( params.bwa_index ){
        bwa_index = file(params.bwa_index)
        if( !bwa_index.exists() ) exit 1, "BWA index not found: $bwa_index"
        bwa_index_amb = file( params.bwa_index+'.amb' )
        if( !bwa_index_amb.exists() ) exit 1, "BWA index not found: $bwa_index_amb"
        bwa_index_ann = file( params.bwa_index+'.ann' )
        if( !bwa_index_ann.exists() ) exit 1, "BWA index not found: $bwa_index_ann"
        bwa_index_bwt = file( params.bwa_index+'.bwt' )
        if( !bwa_index_bwt.exists() ) exit 1, "BWA index not found: $bwa_index_bwt"
        bwa_index_pac = file( params.bwa_index+'.pac' )
        if( !bwa_index_pac.exists() ) exit 1, "BWA index not found: $bwa_index_pac"
        bwa_index_sa = file( params.bwa_index+'.sa' )
        if( !bwa_index_sa.exists() ) exit 1, "BWA index not found: $bwa_index_sa"
    } else if ( params.fasta ){
        fasta = file(params.fasta)
        if( !fasta.exists() ) exit 1, "Fasta file not found: $fasta"
    }
}


params.reads = "data/*{1,2}.fastq.gz"
params.outdir = './results'


// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0



// Header log info
log.info "========================================="
log.info " NouGat : De Novo Best Practice v${version}"
log.info "========================================="


/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }


/*
 * STEP - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    file '*fq.gz' into trimmed_reads_jellyfish, trimmed_reads_bwa
    file '*trimming_report.txt' into trimgalore_results

    script:
    single = reads instanceof Path
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (single) {
        """
        trim_galore --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}

/*
 * STEP - Jellyfish
 */
process jellyfish {
    tag "$reads"
    publishDir "${params.outdir}/jellyfish", mode: 'copy'

    input:
    file reads from trimmed_reads_jellyfish

    output:
    file '*.jf' into jellyfish_db
    file '*.hist' into jellifish_hist

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

    """
    #Getting prefix
    gunzip -c $reads | jellyfish count -o ${prefix}.jf  -m 25 -s 1000M -t ${task.cpus} -C /dev/fd/0
    jellyfish histo -o ${prefix}_jf.hist -f ${prefix}.jf
    """


}



/*
 * STEP - BWA
 */
process bwa {
    tag "$reads"
    publishDir "${params.outdir}/bwa", mode: 'copy'

    when:
    params.align

    input:
    file index from bwa_index
    file bwa_index_sa
    file bwa_index_bwt
    file bwa_index_ann
    file bwa_index_amb
    file bwa_index_pac
    file reads from trimmed_reads_bwa

    output:
    file ('*.bam') into aligned
    stdout into bwa_log

    script:
    """
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
    bwa mem -t ${task.cpus} $index $reads |  samtools sort  --threads $task.cpus - > \$f.bam

    """
}

/*
 * STEP - QUALIMAP
*/
process qualimap {
    tag "$reads"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    params.align

    input:
    file bam from aligned

    output:
    file ('*_stats') into qualimap_result

    script:
    """
    qualimap bamqc -nt 4 -bam $bam

    """
}


/*
 * STEP 5 MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('jellyfish/*') from jellifish_hist.flatten().toList()
    file ('qualimap/*') from qualimap_result.flatten().toList()

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
    which multiqc
    multiqc -f .
    """
}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/


def isAllowedParams(params) {
  final test = true
  params.each{
    if (!checkParams(it.toString().split('=')[0])) {
      println "params ${it.toString().split('=')[0]} is unknown"
      test = false
    }
  }
  return test
}

def checkParams(it) {
  // Check if params is in this given list
  return it in [
    'help',
    'version',
    'genomes',
    'reads',
    'align',
    'bwa_index']
}


def grabRevision() {
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage(version, revision) { // Display help message
  log.info "NouGat de novo pipeline ~ $version - revision: $revision"
  log.info "    usage:"
  log.info "        nextflow run -c configuration.config ...."
  log.info "    --version: prints the version"
  log.info "    --help: prints this help message"
  log.info "    --genomes: still not rsure..."
  log.info "    --reads: ehre the reads are"
  log.info "    --align: if present performs also aligment"
}


def versionMessage(version, revision) { // Display version message
  log.info "NouGat de novo assembly pipeline"
  log.info "  version   : $version"
  log.info workflow.commitId ? "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]" : "  revision  : $revision"
}













