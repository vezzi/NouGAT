import sys, os, yaml, glob
import gzip
import subprocess
import pandas as pd
import re
import string
from de_novo_scilife import align
from de_novo_scilife import common
from de_novo_scilife import QCcontrol

from de_novo_scilife import pdf
from de_novo_scilife.pdf.theme import colors, DefaultTheme


def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if "tools" in sample_config and len(sample_config["tools"]) > 0:
        """If so, execute them one after the other in the specified order (might not work)"""
        for command in sample_config["tools"]:
            """with this I pick up at run time the correct function in the current module"""
            command_fn    = getattr(sys.modules[__name__] , "_run_{}".format(command))
            """Update sample config, each command return sample_config and if necessary it modifies it"""
            sample_config = command_fn(global_config, sample_config, sorted_libraries_by_insert)
    else:
        #run default pipeline for QC
        sample_config = _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_fastqc(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_abyss(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_align(global_config, sample_config, sorted_libraries_by_insert)
    
    #now reverse complement
    #sample_config = _run_revcomp(global_config, sample_config, sorted_libraries_by_insert)
    _run_report(global_config, sample_config,sorted_libraries_by_insert)

def _run_fastqc(global_config, sample_config,sorted_libraries_by_insert):
    return QCcontrol._run_fastqc(global_config, sample_config,sorted_libraries_by_insert)

def _run_abyss(global_config, sample_config,sorted_libraries_by_insert):
    return QCcontrol._run_abyss(global_config, sample_config,sorted_libraries_by_insert)



def _run_align(global_config, sample_config,sorted_libraries_by_insert):
    
    if "reference" not in sample_config:
        print "reference sequence not provided, skypping alignment step. Please provide a reference if you are intrested in aligning the reads against a reference"
        return
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


def _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert):
    program        = global_config["Tools"]["trimmomatic"]["bin"]
    program_folder = os.path.dirname(program)
    if "adapters" not in sample_config:
        sys.exit("running MP pipeline, adapters file to be used in trimming are needed for Trimmomatic. Please specify them\
        in the sample configuration file and rerun")
    adapterFile    = sample_config["adapters"]
    if not os.path.exists(adapterFile):
        sys.exit("Trimmomatic cannot be run as adapter file is not specified or points to unknown position: {}".format(adapterFile))

    mainDirectory   = os.getcwd()
    trimmomaticDir  = os.path.join(mainDirectory, "Trimmomatic")
    if not os.path.exists(trimmomaticDir):
        os.makedirs(trimmomaticDir)
    os.chdir(trimmomaticDir)
    #now I am in running dir, I need to process one by one the libraries
    threads = 8
    if "threads" in sample_config:
        threads = sample_config["threads"]

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if orientation=="outtie": # this must be done only with MP
            if read2 is not None:
                read1_baseName = os.path.split(read1)[1].split(".")[0]
                read2_baseName = os.path.split(read2)[1].split(".")[0]
                output_read1_pair = os.path.join(trimmomaticDir,  "{}.fastq.gz".format(read1_baseName))
                output_read1_sing = os.path.join(trimmomaticDir, "{}_u.fastq.gz".format(read1_baseName))
                output_read2_pair = os.path.join(trimmomaticDir, "{}.fastq.gz".format(read2_baseName))
                output_read2_sing = os.path.join(trimmomaticDir,  "{}_u.fastq.gz".format(read2_baseName))
                command = ["java",  "-jar", program, "PE", "-threads", "{}".format(threads),  "-phred33",  read1, read2,  output_read1_pair ,output_read1_sing , output_read2_pair, output_read2_sing ,"ILLUMINACLIP:{}:2:30:10".format(adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:30" ]
                common.print_command(command)
                if not common.check_dryrun(sample_config) and not os.path.exists("{}.fastq.gz".format(read1_baseName)): # do not execute is files have been already gennerated
                    stdOut = open("{}_trimmomatic.stdOut".format(read1_baseName), "w")
                    stdErr = open("{}_trimmomatic.stdErr".format(read1_baseName), "w")
                    returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr) # run the program
                    if returnValue != 0:
                        print "error while running command: {}".format(command)
                libraryInfo["pair1"]       = output_read1_pair
                libraryInfo["pair2"]       = output_read2_pair
                libraryInfo["trimmomatic"] = os.path.join(trimmomaticDir, "{}_trimmomatic.stdErr".format(read1_baseName))
    os.chdir(mainDirectory)
    return sample_config


def _run_revcomp(global_config, sample_config, sorted_libraries_by_insert):
    currentDir  = os.getcwd()
    workingDir  = os.path.join(currentDir, "reverse_complemented")
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    os.chdir(workingDir)

    common.print_command("reversing complementing reads, internal function")
    for library, libraryInfo in sorted_libraries_by_insert:
        pairs = [libraryInfo["pair1"] , libraryInfo["pair2"]]
        for pair in pairs:
            base_name   = os.path.split(pair)[1].split(".fastq")[0]
            output_name = "{}_rc.fastq".format(base_name)
            if not os.path.exists("{}.gz".format(output_name)):
                rev_compl_file          = _rev_comp_file(pair, output_name)
                libraryInfo["pair1_rc"] = os.path.join(workingDir, rev_compl_file)
    os.chdir(currentDir)
    return sample_config
    
    
def _rev_comp_file(input_fastq, output_fastq):
    proc = subprocess.Popen(["zcat", "{}".format(input_fastq)], stdout=subprocess.PIPE, bufsize=1)
    with gzip.open('{}.gz'.format(output_fastq), 'wb') as output_f:
        for header in iter(proc.stdout.readline, b''):
            header   = header.rstrip()
            sequence = proc.stdout.readline().rstrip()
            comment  = proc.stdout.readline().rstrip()
            quality  = proc.stdout.readline().rstrip()
        
            sequence_rc = rc(sequence)
            quality_rev = reverse(quality)

            output_f.write("{}\n".format(header))
            output_f.write("{}\n".format(sequence_rc))
            output_f.write("{}\n".format(comment))
            output_f.write("{}\n".format(quality_rev))

    proc.communicate()
    return "{}.gz".format(output_fastq)

def reverse(quality):
    seq_rev = quality[::-1]
    return seq_rev


def rc(dna):
    complements = string.maketrans('acgtrymkbdhvnACGTRYMKBDHVN', 'tgcayrkmvhdbnTGCAYRKMVHDBN')
    rcseq = dna.translate(complements)[::-1]
    return rcseq




def _run_report(global_config, sample_config, sorted_libraries_by_insert):
    """This function produces a pdf report and stores the important resutls in a single folder"""
   
    ### retrive all info needed to write the report
    sampleName = "sample"
    if "output" in sample_config:
        sampleName = sample_config["output"]
    projectName = "anonymous_project"
    if "projectName" in sample_config:
        projectName = sample_config["projectName"]


    currentDir  = os.getcwd()
    workingDir  = os.path.join(currentDir, "report")
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    os.chdir(workingDir)

    PDFtitle = "{}.pdf".format(sample_config["output"])

    TABLE_WIDTH = 540 # this you cannot do in rLab which is why I wrote the helper initially
    class MyTheme(DefaultTheme):
        doc = {
            'leftMargin': 25,
            'rightMargin': 25,
            'topMargin': 20,
            'bottomMargin': 25,
            'allowSplitting': False
            }
    # let's create the doc and specify title and author
    doc = pdf.Pdf('{} {}'.format(projectName, sampleName), 'NGI-Stockholm, Science for Life Laboratory')
    
    # now we apply our theme
    doc.set_theme(MyTheme)
    
    # give me some space
    doc.add_spacer()
    
    # this header defaults to H1
    scriptDirectory = os.path.split(os.path.abspath(__file__))[0]
    logo_path = os.path.join(scriptDirectory, '../pictures/SciLifeLab.jpeg')
    doc.add_image(logo_path, 180, 67, pdf.CENTER)
    logo_path = os.path.join(scriptDirectory, '../pictures/NGI.jpeg')
    doc.add_image(logo_path, 180, 67, pdf.CENTER)
    # give me some space
    doc.add_spacer()

    doc.add_header('NGI-Stockholm -- Science For Life Laboratory')
    doc.add_header('Mate Pair Best Practice Analysis Report')
    doc.add_header('{} -- {}'.format(projectName, sampleName))
    # give me some space
    #doc.add_spacer()
    #doc.add_paragraph("NGI-Stockholm and Science For Life Laboratory bla bla bla [Mission statement]")

    doc.add_spacer()
    doc.add_paragraph("For sample {} belonging to project {} NGI-Stockholm best-practice-analysis for Mate Pairs produce with the Nextera MP \
protcol has been run. For more informations about the protocols and the general guidelines for analysis refer to \
http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf".format(sampleName, projectName))
    doc.add_spacer()
    tools = ["trimmomatic", "fastqc", "abyss", "align"]
    if "tools" in sample_config and len(sample_config["tools"]) > 0:
        tools = sample_config["tools"]
    doc.add_paragraph("The following tools have been employed (tools are listed in order of execution):")
    bollet_list = []
    for tool in tools:
        if tool != "align":
            program_path = global_config["Tools"][tool]["bin"]
            bollet_list.append("{} : {}".format(tool, program_path))
        else:
            bollet_list.append("{} : {}".format(tool, global_config["Tools"]["bwa"]["bin"]))
            bollet_list.append("{} : {}".format(tool, global_config["Tools"]["samtools"]["bin"]))
            bollet_list.append("{} : {}".format(tool, global_config["Tools"]["picard"]["bin"]))
    doc.add_list(bollet_list)
    doc.add_spacer()
    doc.add_paragraph("For each tool results are reported (tables, pictures). Moreover you will find all the results and commands that have been\
in the result delivery folder on Uppmax")


    for tool in tools:
        doc.add_header(tool , pdf.H2)
        if tool  == "trimmomatic":
            doc.add_paragraph("As a result of the Illumina Nextera Mate Pair protocol many reads sequenced will contain the adapter sequence. Illumina reccomends\
to remove the adaptor before use the reads in any downstream analysis.")
            doc.add_paragraph("Adapter sequences removed are:")
            adapter_file = sample_config["adapters"]
            adapters     = []
            with open(adapter_file) as file:
                lines       = file.readlines()
                adapters = [lines[1].rstrip(),lines[3].rstrip(),lines[5].rstrip()]
            doc.add_list(adapters)
            doc.add_spacer()

            trimmomatic_table_part1 = [[sampleName, "#orig_pairs", "#survived_pairs"]] # this is the header row
            trimmomatic_table_part2 = [[sampleName,"#survived_fw_only", "#survived_rv_only", "#discarded"]]
            
            total_orig_pairs       = 0
            total_survived_pairs   = 0
            total_survived_fw_only = 0
            total_survived_rv_only = 0
            total_discarded        = 0
            
            for library, libraryInfo in sorted_libraries_by_insert:
                runName = os.path.basename(libraryInfo["trimmomatic"]).split("_1_trimmomatic.stdErr")[0]
                with open(libraryInfo["trimmomatic"]) as trimmomatic_output:
                    lines       = trimmomatic_output.readlines()
                    result_line = lines[-2].rstrip()
                    match_string = re.compile('Input Read Pairs: (\d+) Both Surviving: (\d+) \(.+\) Forward Only Surviving: (\d+) \(.+\) Reverse Only Surviving: (\d+) \(.+\) Dropped: (\d+) \(.+\)')
                    read_pairs       = int(match_string.match(result_line).group(1))
                    survived_pairs   = int(match_string.match(result_line).group(2))
                    survived_fw_only = int(match_string.match(result_line).group(3))
                    survived_rv_only = int(match_string.match(result_line).group(4))
                    discarded        = int(match_string.match(result_line).group(5))
                    read_pairs_perc          =  "({0:.0f}%)".format((float(survived_pairs)/read_pairs) * 100)
                    survived_fw_only_perc    =  "({0:.0f}%)".format((float(survived_fw_only)/read_pairs) * 100)
                    survived_rv_only_perc    =  "({0:.0f}%)".format((float(survived_rv_only)/read_pairs) * 100)
                    survived_discarded_perc  =  "({0:.0f}%)".format((float(discarded)/read_pairs) * 100)

                    total_orig_pairs       += read_pairs
                    total_survived_pairs   += survived_pairs
                    total_survived_fw_only += survived_fw_only
                    total_survived_rv_only += survived_rv_only
                    total_discarded        += discarded
                trimmomatic_table_part1.append([runName,read_pairs, "{} {}".format(survived_pairs, read_pairs_perc)]) # these are the other rows
                trimmomatic_table_part2.append([runName, "{} {}".format(survived_fw_only, survived_fw_only_perc), \
                "{} {}".format(survived_rv_only, survived_rv_only_perc),  "{} {}".format(discarded, survived_discarded_perc)]) # these are the other rows

            survived_pairs_perc            =  "({0:.0f}%)".format((float(total_survived_pairs)/total_orig_pairs) * 100)
            survived_survived_fw_only_perc =  "({0:.0f}%)".format((float(total_survived_fw_only)/total_orig_pairs) * 100)
            survived_survived_rv_only_perc =  "({0:.0f}%)".format((float(total_survived_rv_only)/total_orig_pairs) * 100)
            survived_discarded_perc        =  "({0:.0f}%)".format((float(total_discarded)/total_orig_pairs) * 100)
            
            
            trimmomatic_table_part1.append(["total", total_orig_pairs, "{} {}".format(total_survived_pairs, survived_pairs_perc)]) # last row is the sum
            trimmomatic_table_part2.append(["total", "{} {}".format(survived_fw_only, survived_fw_only_perc), \
                "{} {}".format(survived_rv_only, survived_rv_only_perc),  "{} {}".format(discarded, survived_discarded_perc)]) # last row is the sum
            
            doc.add_table(trimmomatic_table_part1, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(trimmomatic_table_part2, TABLE_WIDTH)
            
            
        if tool == "fastqc":
            fastqc_dir = sample_config["fastqc"]
            for fastqc_run in [dir for dir in os.listdir(fastqc_dir) if os.path.isdir(os.path.join(fastqc_dir, dir))]:
                doc.add_header("{} -- Per Base Quality".format(fastqc_run) , pdf.H3)
                fastqc_run_dir = os.path.join(fastqc_dir, fastqc_run, "Images")
                doc.add_image(os.path.join(fastqc_run_dir,"per_base_quality.png"), 400, 180, pdf.CENTER)
                doc.add_header("{} -- Sequence Length Distribution".format(fastqc_run) , pdf.H3)
                fastqc_run_dir = os.path.join(fastqc_dir, fastqc_run, "Images")
                doc.add_image(os.path.join(fastqc_run_dir,"sequence_length_distribution.png"), 400, 180, pdf.CENTER)
        if tool == "abyss":
            doc.add_paragraph("kmer profile with k={}.".format(sample_config["kmer"]))
            kmer_1_200 = os.path.join(sample_config["abyss"], "kmer_coverage_1_200.png")
            doc.add_image(kmer_1_200, 500, 300, pdf.CENTER)
        if tool == "align":
            alignments       = sample_config["alignments"][0]
            alignment_path   = alignments[1]
            alignment_prefix = alignments[2]
            align_dir        = os.path.split(alignment_path)[0]
            doc.add_header("{} -- Collect Insert Size Metrics".format(sampleName) , pdf.H3)
            with open(os.path.join(align_dir, "{}.collectInsertSize.txt".format(alignment_prefix) )) as collectInsertSize:
                    lines  = collectInsertSize.readlines()
                    line  = lines[6].rstrip().split("\t")
                    insertSize_table = [[line[7], line[6], line[4], line[5]]]  # this is the header row
                    line  = lines[7].rstrip().split("\t")
                    insertSize_table.append([line[7], line[6], line[4], line[5]])  # this is the header row
                    line  = lines[8].rstrip().split("\t")
                    insertSize_table.append([line[7], line[6], line[4], line[5]])  # this is the header row
                    line  = lines[9].rstrip().split("\t")
                    insertSize_table.append([line[7], line[6], line[4], line[5]])  # this is the header row
            doc.add_table(insertSize_table, TABLE_WIDTH)
            doc.add_spacer()
            full_path_to_pdf =  os.path.join(align_dir, "{}.collectInsertSize.pdf".format(alignment_prefix))
            doc.add_paragraph("Insert size plot can be find in the result directory: {}".format(os.path.join(full_path_to_pdf.split(os.sep)[-3], full_path_to_pdf.split(os.sep)[-2], full_path_to_pdf.split(os.sep)[-1])))
            doc.add_spacer()
            
            doc.add_header("{} -- Duplicate Metrics".format(sampleName) , pdf.H3)
            with open(os.path.join(align_dir, "{}.markDuplicates.txt".format(alignment_prefix) )) as collectInsertSize:
                    lines  = collectInsertSize.readlines()
                    line  = lines[6].rstrip().split("\t")
                    duplication_table_part1 = [line[0:3]]  # this is the header row
                    duplication_table_part2 = [line[4:6]]  # this is the header row
                    duplication_table_part3 = [line[7:9]]  # this is the header row
                    line  = lines[7].rstrip().split("\t")
                    duplication_table_part1.append(line[0:3])  #
                    duplication_table_part2.append(line[4:6])  #
                    duplication_table_part3.append(line[7:9])  #
            doc.add_table(duplication_table_part1, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(duplication_table_part2, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(duplication_table_part3, TABLE_WIDTH)
            doc.add_spacer()
            full_path_to_bam =  os.path.join(align_dir, "{}_noDup.bam".format(alignment_prefix))
            doc.add_paragraph("Bam file with makeked duplicated reads can be faound at: {}".format(os.path.join(full_path_to_bam.split(os.sep)[-3], full_path_to_bam.split(os.sep)[-2], full_path_to_bam.split(os.sep)[-1])))
            doc.add_spacer()
            
            
    doc.render(PDFtitle)
    
    os.chdir(currentDir)
    










