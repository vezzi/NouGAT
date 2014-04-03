import sys, os, yaml, glob
import subprocess
import pandas as pd
import gzip
import re
import string
import shutil
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from de_novo_scilife import common
from de_novo_scilife import align


def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    _check_libraries(sorted_libraries_by_insert)
    sample_config = _build_new_reference(sample_config) # filter out short contigs
    if "tools" in sample_config:
        """If so, execute them one after the other in the specified order (might not work)"""
        for command in sample_config["tools"]:
            """with this I pick up at run time the correct function in the current module"""
            command_fn    = getattr(sys.modules[__name__] , "_run_{}".format(command))
            """Update sample config, each command return sample_config and if necessary it modifies it"""
            sample_config = command_fn(global_config, sample_config, sorted_libraries_by_insert)
    else:
        #run default pipeline for QC
        sample_config = _run_align(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_qaTools(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_FRC(global_config, sample_config, sorted_libraries_by_insert)



def _run_align(global_config, sample_config,sorted_libraries_by_insert):
    if "reference" not in sample_config:
        print "reference sequence not provided, skypping alignment step. Please provide a reference if you are intrested in aligning the reads against a reference"
        return sample_config
    if not os.path.exists("alignments"):
        os.makedirs("alignments")
    os.chdir("alignments")
    sorted_libraries_by_insert =  align._align_reads(global_config, sample_config,  sorted_libraries_by_insert) # align reads
    sorted_alignments_by_insert = align._merge_bam_files(global_config, sample_config, sorted_libraries_by_insert) # merge alignments
    sorted_alignments_by_insert = align.picard_CGbias(global_config, sample_config,sorted_alignments_by_insert) # compute picard stats
    sorted_alignments_by_insert = align.picard_collectInsertSizeMetrics(global_config, sample_config,sorted_alignments_by_insert)
    sorted_alignments_by_insert = align.picard_markDuplicates(global_config, sample_config,sorted_alignments_by_insert)
    os.chdir("..")
    sample_config["alignments"] = sorted_alignments_by_insert
    return sample_config


def _check_libraries(sorted_libraries_by_insert):
    different_inserts   = 0
    current_insert      = -1
    orientation         = ""
    for library, libraryInfo in sorted_libraries_by_insert:
        if current_insert == -1:
            current_insert    = libraryInfo["insert"]
            different_inserts = 1
        else :
            if current_insert != libraryInfo["insert"]:
                current_insert    = libraryInfo["insert"]
                different_inserts += 1
        if libraryInfo["orientation"] == "outtie":
            sys.exit("error: in valiadation only innie libraries can be used, please reverse complement your MP libs (N.B. if you employed MP module of this pipeline then the reverse complemented reads are already available). This is needed in order to run FRCurve, please complain with FRCurve author!!!!")
    if different_inserts > 2:
        sys.exit("error: in valiadation only two libraries are admitted (usually a PE and a MP, sometimes 2 PE)")
    return



def _build_new_reference(sample_config):
    minCtgLength = 500
    if "minCtgLength" in sample_config:
        minCtgLength = sample_config["minCtgLength"]
        if minCtgLength < 500:
            sys.exit("min contig length must be higher than 500bp, lower values will complicate the job of valiadation tools and make results difficult to interpret. For mammalian genomes minCtgLength > 1Kbp is strongly suggested")
    reference      = sample_config["reference"]
    reference_dir  = os.path.abspath("reference")
    if not os.path.exists(reference_dir):
        os.makedirs(reference_dir)
    os.chdir(reference_dir)
    new_reference_name = os.path.abspath(os.path.basename(reference))
    if os.path.exists(new_reference_name):
        sample_config["reference"] = new_reference_name
        os.chdir("..")
        return sample_config # already created the new reference
    with open(new_reference_name, "w") as new_ref_fd:
        with open(reference, "r") as ref_fd:
            fasta_header    = ref_fd.readline()
            sequence        = ""
            for line in ref_fd:
                line = line
                if line.startswith(">"):
                    if len(sequence) >= minCtgLength:
                        new_ref_fd.write(fasta_header)
                        new_ref_fd.write(sequence)
                    sequence        = ""
                    fasta_header    = line
                else:
                    sequence+=line
            if len(sequence) >= minCtgLength:
                new_ref_fd.write(fasta_header)
                new_ref_fd.write(sequence)
    sample_config["reference"] = new_reference_name
    os.chdir("..")
    return sample_config


def _run_CEGMA(global_config, sample_config, sorted_alignments_by_insert):
    cegma       = global_config["evaluete"]["CEGMA"]["bin"]
    assembly    = sample_config["reference"]
    #module load cegma/2.4.010312
    return sample_config



def _run_FRC(global_config, sample_config, sorted_libraries_by_insert):
    mainDir       = os.getcwd()
    FRCurveFolder = os.path.join(os.getcwd(), "FRCurve")
    if not os.path.exists(FRCurveFolder):
        os.makedirs(FRCurveFolder)
    os.chdir("FRCurve")
    program=global_config["Tools"]["FRC"]["bin"]

    genomeSize  = sample_config["genomeSize"]
    reference   = sample_config["reference"]
    output      = sample_config["output"]
    alignments  = sample_config["alignments"]
    
    peBam       = alignments[0][1]
    peInsert    = alignments[0][0]
    peMinInsert = int(peInsert - peInsert*0.30)
    peMaxInsert = int(peInsert + peInsert*0.30)
    command = [program, "--pe-sam", peBam, "--pe-min-insert", "{}".format(peMinInsert) , "--pe-max-insert", "{}".format(peMaxInsert), "--CEstats-PE-min", "-5", "--CEstats-PE-max", "5"]
    if len(alignments) > 1:
        mpBam       = alignments[1][1]
        mpInsert    = alignments[1][0]
        mpMinInsert = int(mpInsert - mpInsert*0.40)
        mpMaxInsert = int(mpInsert + mpInsert*0.40)
        command += ["--mp-sam", mpBam, "--mp-min-insert", "{}".format(mpMinInsert), "--mp-max-insert", "{}".format(mpMaxInsert)]
    command += [ "--genome-size", "{}".format(genomeSize), "--output", output]
    common.print_command(command)
    if not common.check_dryrun(sample_config) and not os.path.exists("{}_FRC.png".format(output)):
        stdOut = open("FRC.stdOut", "a")
        stdErr = open("FRC.stdErr", "a")
        returnValue = subprocess.call(command , stdout=stdOut , stderr=stdErr)
        if not returnValue == 0:
            sys.exit("error, while running FRCurve: {}".format(command))
        plotFRCurve(output)
    os.chdir("..")
    return sample_config


def plotFRCurve(output):
    names = ["_FRC" , "COMPR_MP_FRC" , "COMPR_PE_FRC" , "HIGH_COV_PE_FRC" , "HIGH_NORM_COV_PE_FRC" ,"HIGH_OUTIE_MP_FRC" , "HIGH_OUTIE_PE_FRC" , "HIGH_SINGLE_MP_FRC" , "HIGH_SINGLE_PE_FRC" , "HIGH_SPAN_MP_FRC" , "HIGH_SPAN_PE_FRC" ,"LOW_COV_PE_FRC" , "LOW_NORM_COV_PE_FRC" , "STRECH_MP_FRC" , "STRECH_PE_FRC"]
    for name in names:
        FRC_data    = pd.io.parsers.read_csv("{}{}.txt".format(output, name), sep=' ', header=None)
        FRC_features= FRC_data[FRC_data.columns[0]].tolist()
        FRC_coverage= FRC_data[FRC_data.columns[1]].tolist()
        plt.plot(FRC_features, FRC_coverage)
        if name == "_FRC":
            plt.title('Feature Resonse Curve -- All Features')
        else:
            plt.title('Feature Resonse Curve -- {}'.format(name))
        plt.plot(FRC_features, FRC_coverage)
        plt.savefig("{}{}.png".format(output, name))
        plt.clf()
    return 0

def _run_qaTools(global_config, sample_config, sorted_libraries_by_insert):
    mainDir       = os.getcwd()
    qaToolsFolder = os.path.join(os.getcwd(), "QAstats")
    if not os.path.exists(qaToolsFolder):
        os.makedirs(qaToolsFolder)
    os.chdir("QAstats")
    program=global_config["Tools"]["qaTools"]["bin"]

    genomeSize  = sample_config["genomeSize"]
    reference   = sample_config["reference"]
    output      = sample_config["output"]
    alignments  = sample_config["alignments"][0]
    BAMfile     = alignments[1]


    command = ["{}".format(program),  "-m",  "-q", "0", "-i",  BAMfile, "{}.cov".format(os.path.basename(BAMfile))]
    common.print_command(command)
    if not common.check_dryrun(sample_config) and not os.path.exists("{}.cov".format(os.path.basename(BAMfile))):
        stdOut = open("QAtools.stdOut", "a")
        stdErr = open("QAtools.stdErr", "a")
        returnValue = subprocess.call(command , stdout=stdOut , stderr=stdErr)
        if not returnValue == 0:
            sys.exit("error, while running QAtools: {}".format(command))
        #now add GC content
        QAtools_dict = {}
        header       = ""
        with open( "{}.cov".format(os.path.basename(BAMfile)), "r") as QA_csv:
            header = QA_csv.readline().rstrip()
            for line in QA_csv:
                line = line.strip().split("\t")
                QAtools_dict[line[0]] = [line[1],line[2],line[3]]
        QA_GC_file = "{}.cov.gc".format(os.path.basename(BAMfile))
        with open(QA_GC_file, "w") as QA_GC_fd:
            QA_GC_fd.write("{}\tGCperc\n".format(header))
            with open(reference, "r") as ref_fd:
                fasta_raw_header    = ref_fd.readline().strip()
                fasta_raw_header    = fasta_raw_header.split(" ")[0]
                fasta_raw_header    = fasta_raw_header.split("\t")[0]
                fasta_header        = fasta_raw_header.split(">")[1]
                sequence            = ""
                for line in ref_fd:
                    line = line.strip()
                    if line.startswith(">"):
                        GC = computeGC(sequence)
                        if fasta_header not in QAtools_dict:
                            sys.exit("error while parsing QAcompute output: probably some wired contig name is present in your assmebly file")
                        QA_GC_fd.write("{}\t{}\t{}\t{}\t{}\n".format(fasta_header, QAtools_dict[fasta_header][0], QAtools_dict[fasta_header][1], QAtools_dict[fasta_header][2], GC ))
                        sequence = ""
                        fasta_raw_header    = line.split(" ")[0]
                        fasta_raw_header    = fasta_raw_header.split("\t")[0]
                        fasta_header        = fasta_raw_header.split(">")[1]
                    else:
                        sequence+=line
                GC = computeGC(sequence)
                if fasta_header not in QAtools_dict:
                    sys.exit("error while parsing QAcompute output: probably some wired contig name is present in your assmebly file")
                QA_GC_fd.write("{}\t{}\t{}\t{}\t{}\n".format(fasta_header, QAtools_dict[fasta_header][0], QAtools_dict[fasta_header][1], QAtools_dict[fasta_header][2], GC ))
        plotQA(QA_GC_file)
    os.chdir("..")
    return sample_config


def plotQA(QA_GC_file):
    #QA_GC_file="lib_500.bam.cov.gc"
    import shutil as sh
    sh.copy(QA_GC_file, "Contigs_Cov_SeqLen_GC.csv")
    QA_data     = pd.io.parsers.read_csv("Contigs_Cov_SeqLen_GC.csv", sep='\t', header=0)
    GCperc      = QA_data['GCperc'].tolist()
    MedianCov   = QA_data['Median_Cov'].tolist()
    SeqLen      = QA_data['Seq_len'].tolist()
    Mean_MedianCov = sum(MedianCov) / float(len(MedianCov))
    Max_MedianCov  = max(MedianCov)
    if Max_MedianCov > 2.5* Mean_MedianCov:
        Max_MedianCov = Mean_MedianCov*2
    #GC_vs_Median Coverage
    plt.plot(GCperc, MedianCov, 'or')
    plt.title('GC content vs Median Coverage')
    plt.xlabel('%GC')
    plt.ylabel('Coverage')
    plotname = "GC_vs_Coverage.png"
    plt.savefig(plotname)
    plt.clf()
    # GC_vs_median eliminate outliers
    plt.plot(GCperc, MedianCov, 'or')
    plt.ylim((10, Max_MedianCov))
    plt.title('GC content vs Median Coverage')
    plt.xlabel('%GC')
    plt.ylabel('Coverage')
    plotname = "GC_vs_Coverage_noOutliers.png"
    plt.savefig(plotname)
    plt.clf()
    
    #Coverage Distribution Histogram
    n, bins, patches = plt.hist(MedianCov, 100, facecolor='g')
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.title('Coverage Distribution')
    plotname = "Coverage_distribution.png"
    plt.savefig(plotname)
    plt.clf()
    #Coverage Distribution Histogram eliminate outliers
    n, bins, patches = plt.hist(MedianCov, 100, facecolor='g', range=(4,Max_MedianCov))
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.title('Coverage Distribution')
    plotname = "Coverage_distribution_noOutliers.png"
    plt.savefig(plotname)
    plt.clf()
    
    
    #Median Cov vs Sequence Length
    plt.plot(MedianCov, map(lambda x: x/1000, SeqLen), 'ro')
    plt.title('Median Coverage vs Contig Length')
    plt.xlabel('Median Coverage')
    plt.ylabel('Contig Length (Kbp)')
    plotname = "MedianCov_vs_CtgLength.png"
    plt.savefig(plotname)
    plt.clf()
    #Median Cov vs Sequence Length eliminate outliers
    plt.plot(MedianCov, map(lambda x: x/1000, SeqLen), 'ro')
    plt.xlim((10, Max_MedianCov))
    plt.title('Median Coverage vs Contig Length')
    plt.xlabel('Median Coverage')
    plt.ylabel('Contig Length (Kbp)')
    plotname = "MedianCov_vs_CtgLength_noOutliers.png"
    plt.savefig(plotname)
    plt.clf()
    
    
    #GC content vs Contig length
    plt.plot(GCperc, map(lambda x: x/1000, SeqLen), 'ro')
    plt.title('%GC vs Contig Length')
    plt.xlabel('%GC')
    plt.ylabel('Contig Length (Kbp)')
    plotname = "GC_vs_CtgLength.png"
    plt.savefig(plotname)
    plt.clf()
    return 0


def computeGC(sequence):
    gcCount         = len(re.findall("[GC]", sequence)) + len(re.findall("[gc]", sequence))
    totalBaseCount  = len(re.findall("[GCTA]", sequence)) + len(re.findall("[gcta]", sequence))
    gcFraction = float(gcCount) / totalBaseCount
    return gcFraction

