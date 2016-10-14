from __future__ import absolute_import
from __future__ import print_function
import sys, os, yaml, glob
from yaml import YAMLError
import subprocess
import argparse
import itertools
import re
import shutil
from nougat import common, align, pdf
from nougat.pdf.theme import colors, DefaultTheme
from nougat.pdf.peakdetect import peakdet


def main(args):
    if "delivery_folder" not in args:
        delivery_folder = os.getcwd() #stage place is the local directory
    delivery_folder = os.path.abspath(args.delivery_folder)
    for qdir in os.listdir(args.qc_folder):
        # Dumping the pipeline state to yaml.. Not elegant, but gets the job done
        sample_yaml = os.path.join(args.qc_folder, qdir, "{}.nougat".format(qdir))
        sample_yaml = os.path.realpath(sample_yaml)
        global_yaml = os.path.realpath(args.global_config)
        try:
            os.chdir(os.path.join(args.qc_folder, qdir))
            with open(sample_yaml) as sample_config_handle:
                sample_config = yaml.load(sample_config_handle)
            with open(global_yaml) as global_config_handle:
                global_config = yaml.load(global_config_handle)
        except IOError as e:
            print("Cannot open file: {}".format(e))
        except YAMLError as e:
            print("Error in config file: {}".format(e))
        else:
            _run_qc_report(global_config, sample_config, delivery_folder)


def _run_qc_report(global_config, sample_config, delivery_folder):
    """This function produces a pdf report and stores the important \
            resutls in a single folder"""

    sorted_libraries_by_insert = common._sort_libraries_by_insert(
            sample_config)
    ### retrive all info needed to write the report
    sampleName = "sample"
    if "output" in sample_config:
        sampleName = sample_config["output"]
    projectName = "anonymous_project"
    if "projectName" in sample_config:
        projectName = sample_config["projectName"]
    
    currentDir = os.getcwd()
    workingDir = os.path.join(currentDir, sampleName)
    #create delivery dir for this sample
    sample_delivery_dir = os.path.join(delivery_folder, sampleName)
    if not os.path.exists(sample_delivery_dir):
        os.makedirs(sample_delivery_dir)

    reportDir   = os.path.join(sample_delivery_dir, "report")
    if not os.path.exists(reportDir):
        os.makedirs(reportDir)

    PDFtitle = os.path.join(sample_delivery_dir, "report",
            "{}.pdf".format(sample_config["output"]))

    # this you cannot do in rLab which is why I wrote the helper initially
    TABLE_WIDTH = 540
    class MyTheme(DefaultTheme):
        doc = {
            'leftMargin': 25,
            'rightMargin': 25,
            'topMargin': 20,
            'bottomMargin': 25,
            'allowSplitting': False
            }
    # let's create the doc and specify title and author
    doc = pdf.Pdf('{} {}'.format(projectName, sampleName),
            'NGI-Stockholm, Science for Life Laboratory')

    # now we apply our theme
    doc.set_theme(MyTheme)
    # give me some space
    doc.add_spacer()
    # this header defaults to H1
    scriptDirectory = os.path.split(os.path.abspath(__file__))[0]
    logo_path = os.path.join(scriptDirectory, '../pictures/ngi_scilife.png')
    doc.add_image(logo_path, 540, 50, pdf.CENTER)
    # give me some space
    doc.add_spacer()

    doc.add_header('NGI-Stockholm -- Science For Life Laboratory')
    doc.add_header('Best-practice analysis for quality checking report')
    doc.add_header('{} -- {}'.format(projectName, sampleName))
    # give me some space
    doc.add_spacer()
    doc.add_paragraph("For sample {} belonging to the project {} "
            "NGI-Stockholm best-practice analysis for quality checking has "
            "been performed. For mate pair libraries produced with Nextera, "
            "best-practice analysis described at this address has been "
            "performed: http://res.illumina.com/documents/products/technotes/"
            "technote_nextera_matepair_data_processing.pdf".format(sampleName,
            projectName))
    doc.add_spacer()
    tools = ["trimmomatic", "fastqc", "abyss", "align", "kmergenie"]
    if "tools" in sample_config and len(sample_config["tools"]) > 0:
        tools = sample_config["tools"]
    doc.add_paragraph("The following tools have been employed \
            (tools are listed in order of execution):")
    bollet_list = []
    for tool in tools:
        if tool != "align":
            program_path = global_config["Tools"][tool]["bin"]
            bollet_list.append("{} : {}".format(tool, program_path))
        else:
            bollet_list.append("{} : {}".format(tool, \
                    global_config["Tools"]["bwa"]["bin"]))
            bollet_list.append("{} : {}".format(tool, \
                    global_config["Tools"]["samtools"]["bin"]))
            bollet_list.append("{} : {}".format(tool, \
                    global_config["Tools"]["picard"]["bin"]))
    doc.add_list(bollet_list)
    doc.add_spacer()
    doc.add_paragraph("The results from each tool is reported in the "
            "following sections. Moreover you will find all the results and "
            "commands that have been run in the delivery folder on Uppmax")

    for tool in tools:
        doc.add_pagebreak()
        doc.add_header(tool.title() , pdf.H2)
        if tool  == "trimmomatic":
            doc.add_paragraph("Reads (both paired and mate pairs) can "
                    "contain parts of the adapter sequence or, in the case of "
                    "mate pairs, part of the linker sequence. Illumina "
                    "recommends to remove the adapter before use of the reads "
                    "in any downstream analysis (this is mandatory for mate "
                    "pairs).")
            doc.add_paragraph("Adapter sequences removed are:")
            adapter_file = sample_config["adapters"]
            adapters     = []
            with open(adapter_file) as afile:
                lines       = afile.readlines()
                for index in range(1, len(lines), 2):
                    adapters.append(lines[index].rstrip())
            doc.add_list(adapters)
            doc.add_spacer()

            trimmomatic_table_part1 = [[sampleName, "#orig_pairs",
                    "#survived_pairs"]] # this is the header row
            trimmomatic_table_part2 = [[sampleName,"#survived_fw_only",
                    "#survived_rv_only", "#discarded"]]

            total_orig_pairs = 0
            total_survived_pairs = 0
            total_survived_fw_only = 0
            total_survived_rv_only = 0
            total_discarded = 0

            for library, libraryInfo in sorted_libraries_by_insert:
                runName = os.path.basename(libraryInfo["trimmomatic"]).split(
                        "_1_trimmomatic.stdErr")[0]
                with open(libraryInfo["trimmomatic"]) as trimmomatic_output:
                    lines       = trimmomatic_output.readlines()
                    result_line = lines[-2].rstrip()
                    match_string = re.compile("Input Read Pairs: (\d+) Both "
                            "Surviving: (\d+) \(.+\) Forward Only Surviving: "
                            "(\d+) \(.+\) Reverse Only Surviving: (\d+) \(.+\) "
                            "Dropped: (\d+) \(.+\)")
                    read_pairs = int(match_string.match(result_line).group(1))
                    survived_pairs = int(match_string.match(
                        result_line).group(2))
                    survived_fw_only = int(match_string.match(
                        result_line).group(3))
                    survived_rv_only = int(match_string.match(
                        result_line).group(4))
                    discarded        = int(match_string.match(
                        result_line).group(5))
                    read_pairs_perc =  "({0:.0f}%)".format(
                            (float(survived_pairs)/read_pairs) * 100)
                    survived_fw_only_perc = "({0:.0f}%)".format(
                            (float(survived_fw_only)/read_pairs) * 100)
                    survived_rv_only_perc = "({0:.0f}%)".format(
                            (float(survived_rv_only)/read_pairs) * 100)
                    survived_discarded_perc = "({0:.0f}%)".format(
                            (float(discarded)/read_pairs) * 100)

                    total_orig_pairs += read_pairs
                    total_survived_pairs += survived_pairs
                    total_survived_fw_only += survived_fw_only
                    total_survived_rv_only += survived_rv_only
                    total_discarded += discarded
                # these are the other rows
                trimmomatic_table_part1.append([runName,read_pairs,
                    "{} {}".format(survived_pairs, read_pairs_perc)])
                trimmomatic_table_part2.append([runName,
                    "{} {}".format(survived_fw_only, survived_fw_only_perc),
                    "{} {}".format(survived_rv_only, survived_rv_only_perc),
                    "{} {}".format(discarded, survived_discarded_perc)])
            survived_pairs_perc =  "({0:.0f}%)".format(
                    (float(total_survived_pairs)/total_orig_pairs) * 100)
            survived_survived_fw_only_perc =  "({0:.0f}%)".format(
                    (float(total_survived_fw_only)/total_orig_pairs) * 100)
            survived_survived_rv_only_perc =  "({0:.0f}%)".format(
                    (float(total_survived_rv_only)/total_orig_pairs) * 100)
            survived_discarded_perc =  "({0:.0f}%)".format(
                    (float(total_discarded)/total_orig_pairs) * 100)
            trimmomatic_table_part1.append(["total", total_orig_pairs,
                "{} {}".format(total_survived_pairs, survived_pairs_perc)])
            # last row is the sum
            trimmomatic_table_part2.append(["total", "{} {}".format(
                survived_fw_only, survived_fw_only_perc),
                "{} {}".format(survived_rv_only, survived_rv_only_perc),
                "{} {}".format(discarded, survived_discarded_perc)])
            doc.add_table(trimmomatic_table_part1, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(trimmomatic_table_part2, TABLE_WIDTH)
            ##now save the trimmed reads
            trimmomaticDir = os.path.split(libraryInfo["trimmomatic"])[0]
            trimmomaticResultDir  = os.path.join(sample_delivery_dir, "fastq_trimmed")
            if not os.path.exists(trimmomaticResultDir):
                os.makedirs(trimmomaticResultDir)
            filesToCopy = [os.path.join(trimmomaticDir, f) for f in \
                    os.listdir(trimmomaticDir) \
                    if (os.path.isfile(os.path.join(trimmomaticDir,f)) \
                    and re.search('.gz$',f))]
            for source in filesToCopy:
                dest = os.path.join(trimmomaticResultDir , os.path.split(source)[1])
                if not os.path.isfile(dest):
                    shutil.copyfile(source, dest)

        if tool == "fastqc" and "fastqc" in sample_config:
            fastqc_dir = sample_config["fastqc"]
            for fastqc_run in [dir for dir in os.listdir(fastqc_dir) \
                    if os.path.isdir(os.path.join(fastqc_dir, dir))]:
                fastqc_run_dir = os.path.join(fastqc_dir, fastqc_run, "Images")
                doc.add_image(os.path.join(fastqc_run_dir,
                    "per_base_quality.png"), 400, 180, pdf.CENTER,
                    "{} -- Per Base Quality".format(fastqc_run))
                fastqc_run_dir = os.path.join(fastqc_dir, fastqc_run, "Images")
                doc.add_image(os.path.join(fastqc_run_dir,
                    "sequence_length_distribution.png"), 400, 180, pdf.CENTER,
                    "{} -- Sequence Length Distribution".format(fastqc_run))
            #If I have not yet copied fastqc results do it
            fastqcResultDir  = os.path.join(sample_delivery_dir, "fastqc")
            if not os.path.exists(fastqcResultDir):
                os.makedirs(fastqcResultDir)
            dirsToBeCopied = [os.path.join(fastqc_dir, f) for f in \
                    os.listdir(fastqc_dir) \
                    if os.path.isdir(os.path.join(fastqc_dir, f))]
            for source in dirsToBeCopied:
                dest = os.path.join(fastqcResultDir, os.path.split(source)[1])
                if not os.path.exists(dest):
                    shutil.copytree(source, dest)

        if tool == "abyss" and "abyss" in sample_config:
            doc.add_paragraph("A possible way to assess the complexity of a "
                    "library even in absence of a reference sequence is to "
                    "look at the kmer profile of the reads. The idea is to "
                    "count all the kmers (i.e., sequence of length k) that occur "
                    "in the reads. In this way it is possible to know how many "
                    "kmers occur 1,2,..., N times and represent this as a "
                    "plot. This plot tell us for each x, how many k-mers "
                    "(y-axis) are present in the dataset in exactly x-copies. "
                    "In an ideal world (no errors in sequencing, no bias, no "
                    "repeating regions) this plot should be as close as "
                    "possible to a gaussian distribution. In reality we will "
                    "always see a peak for x=1 (i.e., the errors) and another "
                    "peak close to the expected coverage. If the genome is "
                    "highly heterozygous a second peak at half of the coverage "
                    "can be expected.")
            kmer_1_200 = os.path.join(sample_config["abyss"],
                    "kmer_coverage.png")
            doc.add_image(kmer_1_200, 500, 300, pdf.CENTER,
                    "kmer profile with k={}.".format(sample_config["kmer"]))
            #copy the results in resutls
            if not os.path.exists("kmer_analysis"):
                os.mkdir("kmer_analysis")
            kmerDir = sample_config["abyss"]
            filesToCopy = [os.path.join(kmerDir, f) for f in \
                    os.listdir(kmerDir) \
                    if (os.path.isfile(os.path.join(kmerDir,f)) \
                    and re.search('.png$',f))]
            filesToCopy.append(os.path.join(kmerDir, "histogram.hist"))
            abyssResultDir = os.path.join(sample_delivery_dir, "kmer_analysis")
            if not os.path.exists(abyssResultDir):
                os.makedirs(abyssResultDir)
            for source in filesToCopy:
                dest = os.path.join(abyssResultDir, os.path.split(source)[1])
                if not os.path.exists(dest):
                    shutil.copyfile(source, dest)

        if tool == "align" and "alignments" in sample_config:
            alignments = sample_config["alignments"][0]
            alignment_path = alignments[1]
            alignment_prefix = alignments[2]
            align_dir = os.path.split(alignment_path)[0]
            doc.add_header("{} -- Collect Insert Size Metrics".format(
                sampleName) , pdf.H3)
            with open(os.path.join(align_dir,
                "{}.collectInsertSize.txt".format(alignment_prefix))) \
                as collectInsertSize:
                    lines  = collectInsertSize.readlines()
                    line  = lines[6].rstrip().split("\t")
                    # this is the header row
                    insertSize_table = [[line[7], line[6], line[4], line[5]]]
                    line  = lines[7].rstrip().split("\t")
                     # this is the header row
                    insertSize_table.append([line[7], line[6], line[4],
                        line[5]])
                    line  = lines[8].rstrip().split("\t")
                     # this is the header row
                    insertSize_table.append([line[7], line[6], line[4],
                        line[5]])
                    line  = lines[9].rstrip().split("\t")
                     # this is the header row
                    insertSize_table.append([line[7], line[6], line[4],
                        line[5]])
                    doc.add_table(insertSize_table, TABLE_WIDTH)
            doc.add_spacer()
            full_path_to_pdf =  os.path.join(align_dir,
                    "{}.collectInsertSize.pdf".format(alignment_prefix))
            doc.add_paragraph("Insert size plot can be found in the result \
                    directory: {}".format(os.path.join("alignments",
                    "{}.collectInsertSize.pdf".format(alignment_prefix))))
            doc.add_spacer()
            doc.add_header("{} -- Duplicate Metrics".format(sampleName),
                    pdf.H3)
            with open(os.path.join(align_dir,
                "{}.markDuplicates.txt".format(alignment_prefix))) as \
                collectInsertSize:
                    lines  = collectInsertSize.readlines()
                    line  = lines[6].rstrip().split("\t")
                    # this is the header row
                    duplication_table_part1 = [line[0:3]]
                    duplication_table_part2 = [line[4:6]]
                    duplication_table_part3 = [line[7:9]]
                    line  = lines[7].rstrip().split("\t")
                    duplication_table_part1.append(line[0:3])
                    duplication_table_part2.append(line[4:6])
                    duplication_table_part3.append(line[7:9])
            doc.add_table(duplication_table_part1, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(duplication_table_part2, TABLE_WIDTH)
            doc.add_spacer()
            doc.add_table(duplication_table_part3, TABLE_WIDTH)
            doc.add_spacer()
            full_path_to_bam =  os.path.join(align_dir,
                    "{}_noDup.bam".format(alignment_prefix))
            doc.add_paragraph("Bam file with marked duplicate reads can be \
                    found at: {}".format(os.path.join("alignments",
                    "{}_noDup.bam".format(alignment_prefix))))
            doc.add_spacer()
            #copy the results in resutls
            if not os.path.exists("alignments"):
                os.mkdir("alignments")
            filesToCopy = [os.path.join(align_dir, f) for f in \
                    os.listdir(align_dir) \
                    if (os.path.isfile(os.path.join(align_dir,f)) \
                    and re.search('{}'.format(alignment_prefix),f))]
            alignmentResultDir = os.path.join(sample_delivery_dir, "alignments")
            if not os.path.exists(alignmentResultDir):
                os.makedirs(alignmentResultDir)
            for source in filesToCopy:
                dest = os.path.join(alignmentResultDir, os.path.split(source)[1])
                if not os.path.exists(dest):
                    shutil.copyfile(source, dest)

        if tool == "kmergenie" and "kmergenie" in sample_config:
            doc.add_paragraph("Assemblers using a de Bruijn graph strategy "
                    "for contig construction (such as Velvet, ABySS and "
                    "SOAPdenovo) fractures the reads into k-sized substrings "
                    "(k-mers). The k-mer size is vital for the performance of "
                    "these assemblers, and is usually selected considering "
                    "several trade-offs between the size and accuracy of the "
                    "produced contigs. Some assemblers choose the k-mer size "
                    "automatically or builds several assemblies (using "
                    "different k-mers) and / or relies on user input. "
                    "Kmergenie is a lightweight program that suggests a best "
                    "k-mer size based on their relative abundance in the "
                    "genomic reads.")
            kmerdir = sample_config["kmergenie"]
            doc.add_image(os.path.join(kmerdir,"histograms.dat.png"), 400, 300, 
                    pdf.CENTER, ("The plot should be roughly concave and have "
                            "a clear global maximum, if not the predicted best "
                            "k is likely to be inaccurate"))
            #copy everything to results
            kmergenieResultDir =  os.path.join(sample_delivery_dir, "kmergenie")
            dest = kmergenieResultDir
            if not os.path.exists(dest):
                shutil.copytree(kmerdir, dest)
    doc.render(PDFtitle)
    
    # Copy the pipeline files and commands run to the report directory
    filesToCopy = glob.glob(currentDir+"/{}_QCcontrol.*".format(sampleName))
    for cfile in filesToCopy:
        shutil.copyfile(cfile, os.path.join(reportDir, os.path.basename(cfile)))

    with open(os.path.join(reportDir, "commands.txt"), "w") as f:
        f.write(sample_config.get("commands", ""))

    os.chdir(currentDir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--global-config', type=str, required=True,
            help="Global configuration file")
    parser.add_argument('--qc-folder', type=str, required=True,
            help=("Path to where denovo QC pipeline was run. Where you find, "
                "the sample folders P????_??? eg. in, /foo/B.Baz_15_01/01-QC/"))
    parser.add_argument('--delivery-folder', type=str, required=False,
            help="Where the delivery will be staged")
    args = parser.parse_args()
    main(args)

