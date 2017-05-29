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

/*
 * SET UP CONFIGURATION VARIABLES
*/


// Pipeline version
version = 0.1

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
    if( !params.bwa_index && params.align ){
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


params.reads = "*{R1,R2}*.fastq.gz"
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
 * STEP 1 - FastQC
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
 * STEP 2 - Trim Galore!
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
 * STEP 3 - Jellyfish
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
    """
    #Getting prefix
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
    zcat $reads | jellyfish count -o \$f.jf  -m 25 -s 1000M -t ${task.cpus} -C /dev/fd/0
    jellyfish histo -o \$f.hist -f \$f.jf
    """


}



/*
 * STEP 4 BWA
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
    bwa mem -t ${task.cpus} $index $reads | samtools view -Sb - | samtools sort - \$f

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
    file ('jellyfish/*hist') from jellifish_hist.flatten().toList()


    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
    multiqc -f .
    """
}
























