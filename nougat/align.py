import sys, os, yaml, glob
import subprocess
import string
import sys
import gzip
import pandas as pd
import common


def _align_reads(global_config, sample_config, sorted_libraries_by_insert):
    reference =  build_reference_bwa(global_config, sample_config)
    #now prepare for alignment
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        threads     = 8
        dryrun      = common.check_dryrun(sample_config)
        if "threads" in sample_config:
            threads  = sample_config["threads"]
        libraryInfo["alignment"] = align_bwa_mem(global_config, read1, read2,
                reference, threads, dryrun)
    return sorted_libraries_by_insert


def _merge_bam_files(global_config, sample_config, sorted_libraries_by_insert):
    BAMfiles = {};
    reference = sample_config["reference"]

    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the "
                "path and not in global config, please make sure to install "
                "bwa properly")

    numInserts = 0
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        alignment = libraryInfo["alignment"]
        if insert not in BAMfiles:
            BAMfiles[insert] = [alignment]
            numInserts += 1
        else:
            BAMfiles[insert].append(alignment)

    BAMfilesMerged = {}
    for insert, insertGroup in BAMfiles.iteritems():
        dir_insert = "lib_{}".format(insert)
        if numInserts == 1:
            dir_insert = sample_config["output"]
        if not os.path.exists(dir_insert):
            os.makedirs(dir_insert)
        os.chdir(dir_insert)
        #check if file is already present
        bamMerged = "lib_{}.bam".format(insert)
        if numInserts == 1:
            bamMerged = "{}.bam".format(sample_config["output"])

        if os.path.exists(bamMerged):
            BAMfilesMerged[insert] = [os.path.abspath(bamMerged), dir_insert]
            os.chdir("..")
            continue # nothiing to be done for this insert

        if len(insertGroup) == 1: # only one sample file for this insert length
            cl = ["ln", "-s", insertGroup[0], bamMerged]
            returnValue = subprocess.call(cl)
            if  not returnValue == 0:
                sys.exit("error, while soft linking {}".format(insertGroup[0]))
        else:
            command = [samtools, "merge",bamMerged]
            for bamfile in insertGroup:
                command.append(bamfile)

            common.print_command(command)
            returnValue = 0
            if not common.check_dryrun(sample_config):
                returnValue = subprocess.call(command)
                if  not returnValue == 0:
                    sys.exit("error, while merging files {}".format(
                        insertGroup))
        BAMfilesMerged[insert] = [os.path.abspath(bamMerged), dir_insert]
        os.chdir("..")

    sorted_alignments_by_insert = []
    for key in sorted(BAMfilesMerged.iterkeys()):
        sorted_alignments_by_insert.append([key, BAMfilesMerged[key][0],
            BAMfilesMerged[key][1]]) # memorise insert length, bam file, folder
    return sorted_alignments_by_insert



def picard_CGbias(global_config, sample_config, sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile, working_dir in sorted_alignments_by_insert:
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        command= ["java", "-Xmx16g", "-XX:PermSize=2g", "-jar",
                os.path.join(picard, "CollectGcBiasMetrics.jar"),
                "REFERENCE_SEQUENCE={}".format(sample_config["reference"]),
                "INPUT={}".format(BAMfile), \
                "OUTPUT={}.collectGcBias.txt".format(output_header),
                "CHART_OUTPUT={}.collectGcBias.pdf".format(output_header),
                "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT",
                "TMP_DIR=$TMPDIR"]
        returnValue = 0;
        common.print_command(command)
        if not os.path.exists("{}.collectGcBias.pdf".format(output_header)):
            if not common.check_dryrun(sample_config):
                stdOut = open("collectGcBias.stdOut", "w")
                stdErr = open("collectGcBias.stdErr", "w")
                returnValue = subprocess.call(command, stdout=stdOut,
                        stderr=stdErr)
                if not returnValue == 0:
                    print "problem running collectGCBias"
        os.chdir("..")
    return sorted_alignments_by_insert



def picard_collectInsertSizeMetrics(global_config, sample_config,
        sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile, working_dir in sorted_alignments_by_insert:
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        histWide = library * 2
        command= ["java", "-Xmx16g", "-XX:PermSize=2g", "-jar",
                os.path.join(picard, "CollectInsertSizeMetrics.jar"),
                "INPUT={}".format(BAMfile), "MINIMUM_PCT=0",
                "HISTOGRAM_FILE={}.collectInsertSize.pdf".format(
                output_header),
                "OUTPUT={}.collectInsertSize.txt".format(output_header),
                "HISTOGRAM_WIDTH={}".format(histWide),
                "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=$TMPDIR"]
        returnValue = 0;
        common.print_command(command)
        if not os.path.exists("{}.collectInsertSize.pdf".format(
            output_header)):
            if not common.check_dryrun(sample_config):
                stdOut = open("collectInsertSize.stdOut", "w")
                stdErr = open("collectInsertSize.stdErr", "w")
                returnValue = subprocess.call(command, stdout=stdOut,
                        stderr=stdErr)
                if not returnValue == 0:
                    print "problem running CollectInsertSizeMetrics"
        os.chdir("..")
    return sorted_alignments_by_insert


def picard_markDuplicates(global_config, sample_config,
        sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile, working_dir in sorted_alignments_by_insert:
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        command= ["java", "-Xmx16g", "-XX:PermSize=3g", "-jar",
                os.path.join(picard, "MarkDuplicates.jar"),
                "INPUT={}".format(BAMfile), "OUTPUT={}_noDup.bam".format(
                output_header),"METRICS_FILE={0}.markDuplicates.txt".format(
                output_header), "ASSUME_SORTED=true",
                "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=$TMPDIR"]
        returnValue = 0;
        common.print_command(command)
        if not os.path.exists("{}.markDuplicates.txt".format(output_header)):
            if not common.check_dryrun(sample_config):
                stdOut = open("removeDup.stdOut", "w")
                stdErr = open("removeDup.stdErr", "w")
                returnValue = subprocess.call(command, stdout=stdOut,
                        stderr=stdErr)
                if not returnValue == 0:
                    print "problem running MarkDuplicates"
        os.chdir("..")
    return sorted_alignments_by_insert


def build_reference_bwa(global_config, sample_config):
    #build the reference if not available
    reference   = sample_config["reference"]
    program = "bwa"
    if "bwa" in global_config["Tools"]:
        program = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa index: bwa not present in "
                "the path and not in global config, please make sure to "
                "install bwa properly")
    # check if reference provided exisists
    reference = os.path.abspath(reference)
    path_name, base_name  = os.path.split(reference)
    index_path = os.path.join(base_name, "bwa", "{}.bwt".format(reference))
    # check if I have already the bwt index
    if os.path.exists(index_path):
        #index already present, nothing to do
        return reference
    #otherwise I need to build the reference, in this case I build it locally
    if not os.path.exists(reference):
        sys.exit("error, reference file {} does not exists".format(reference))
    # check if bwa index already created
    current_dir           = os.getcwd() 
    bwa_index_folder      = os.path.join(path_name, "bwa")
    #if needed create directory
    if not os.path.exists(bwa_index_folder):
        os.makedirs(bwa_index_folder)
    os.chdir(bwa_index_folder)
    # if needed soft link the reference
    if not os.path.exists(base_name):
        #check and remove broken links
        if os.path.lexists(base_name):
            os.remove(base_name)
        returnValue = subprocess.call(["ln", "-s", reference, base_name])
        if not returnValue == 0:
            sys.exit("error while trying to soft link reference sequence")
    # now I have a soflinked copy
    reference = os.path.join(path_name, "bwa", base_name)
    # now check if index alredy build or not
    if not os.path.exists("{}.bwt".format(reference)):
        # then create the index sequence
        bwa_stdOut = open("bwa_index.stdOut", "w")
        bwa_stdErr = open("bwa_index.stdErr", "w")
        command = [program, "index", reference]
        common.print_command(command)
        if not common.check_dryrun(sample_config):
            returnValue = subprocess.call(command, stdout=bwa_stdOut,
                    stderr=bwa_stdErr)
            if  not returnValue == 0:
                sys.exit("error, while indexing reference file {} "
                        "with bwa index".format(reference))
    #extra control to avoid problem with unexpected return value
    if not os.path.exists("{}.bwt".format(reference)):
        sys.exit("bwa index failed")
    os.chdir(current_dir)
    return reference




def align_bwa_mem(global_config, read1, read2, reference, threads, dryrun):
    aligner = "bwa"
    if "bwa" in global_config["Tools"]:
        aligner = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa mem: bwa not present in the "
                "path and not in global config, please make sure to install "
                "bwa properly")

    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the "
                "path and not in global config, please make sure to install "
                "bwa properly")

    # extract base name
    libraryBase = ""
    if read2:
        libraryBase = os.path.basename(read1).split("_1.fastq")[0]
    else:
        libraryBase = os.path.basename(read1).split(".fastq")[0]

    if not os.path.exists(libraryBase):
        os.makedirs(libraryBase)
    os.chdir(libraryBase)
    mappingBase = "{}_to_{}".format(libraryBase,
            os.path.basename(reference).split(".fasta")[0])
    BAMsorted   = "{}.bam".format(mappingBase)
    BAMunsorted = "{}.unsorted.bam".format(mappingBase)
    SAMMapped   = "{}.unsorted.sam".format(mappingBase)
    if os.path.exists(os.path.abspath(BAMsorted)):
        BAMsorted = os.path.abspath(BAMsorted)
        os.chdir("..")
        return BAMsorted


    bwa_mem_command = [aligner, "mem", "-M", "-t", "{}".format(threads),
            reference, read1, read2]
    samtools_view_command = [samtools, "view", "-b", "-S",  "-u",  "-"]

    if not os.path.exists(BAMunsorted):
        command = "{} | {} > {}".format(" ".join(bwa_mem_command),
                " ".join(samtools_view_command), BAMunsorted)
        bwa_stdOut       = open("bwa.stdOut", "w")
        bwa_stdErr       = open("bwa.stdErr", "w")
        common.print_command(command)
        if not dryrun:
            subprocess.call(command, shell=True, stdout=bwa_stdOut,
                    stderr=bwa_stdErr)

    samtools_sort_command = [samtools, "sort", "-@", "{}".format(threads),
            "-m" , "1G", BAMunsorted,  mappingBase]
    command = " ".join(samtools_sort_command)
    if not os.path.exists(BAMsorted):
        stdOut       = open("sam_sort.stdOut", "w")
        stdErr       = open("sam_sort.stdErr", "w")
        common.print_command(command)
        if not dryrun:
            subprocess.call(command, shell=True, stdout=stdOut, stderr=stdErr)

    if os.path.exists(BAMsorted) and os.path.exists(BAMunsorted):
        subprocess.call(["rm", BAMunsorted])
    BAMsorted = os.path.abspath(BAMsorted)
    os.chdir("..")
    return BAMsorted


def align_bwa(read1, read2, reference):
    aligner = "bwa"
    if "bwa" in global_config["Tools"]:
        program = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa aln: bwa not present in the "
                "path and not in global config, please make sure to install "
                "bwa properly")

    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        program = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the "
                "path and not in global config, please make sure to install "
                "bwa properly")

    libraryBase = os.path.basename(read1).split(".")[0]
    if read2:
        libraryBase = "{}_{}".format(libraryBase.split("_")[0],
                libraryBase.split("_")[1])
    if not os.path.exists(libraryBase):
        os.makedirs(libraryBase)
    os.chdir(libraryBase)
    mappingBase = "{}_to_{}".format(libraryBase,
            os.path.basename(reference).split(".fasta")[0])
    BAMsorted   = "{}.bam".format(mappingBase)
    BAMunsorted = "{}.unsorted.bam".format(mappingBase)
    if os.path.exists(os.path.abspath(BAMsorted)):
        BAMsorted = os.path.abspath(BAMsorted)
        os.chdir("..")
        return BAMsorted

    if not os.path.exists(BAMunsorted):
        bwaAln_1 = "{}_1.sai".format(mappingBase);
        bwaAln_2 = "{}_2.sai".format(mappingBase);
        if not os.path.exists("{}_1.sai".format(mappingBase)):
            with open(bwaAln_1, "w") as fh:
                returnValue = subprocess.call(["bwa", "aln", "-t", "8",
                    reference, read1], stdout=fh)
                if  not returnValue == 0:
                    sys.exit("error, while aligning read {} \
                            against {}".format(read1, reference))
        if read2:
            if not os.path.exists("{}_2.sai".format(mappingBase)):
                with open(bwaAln_2, "w") as fh:
                    returnValue = subprocess.call(["bwa", "aln", "-t",
                        "8", reference, read2], stdout=fh)
                    if  not returnValue == 0:
                        sys.exit("error, while aligning read {} \
                                against {}".format(read1, reference))
        with open(BAMunsorted, 'w') as fh:
            if read2:
                cl1 = ["bwa", "sampe", "-P", "-s",  reference, bwaAln_1,
                        bwaAln_2, read1, read2]
                cl2 = ["samtools", "view", "-Shb", "-F", "4", "-"]
                p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cl2, stdin=p1.stdout,
                        stdout=fh).communicate()
                p1.stdout.close()
            else:
                cl1 = ["bwa", "samse",  reference, bwaAln_1, read1]
                cl2 = ["samtools", "view", "-Shb", "-F", "4", "-"]
                p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cl2, stdin=p1.stdout,
                        stdout=fh).communicate()
                p1.stdout.close()

    subprocess.call(["rm", bwaAln_1, bwaAln_2])
    BAMsortedHeader = "{}".format(mappingBase)
    subprocess.call(["samtools", "sort",  BAMunsorted, BAMsortedHeader])
    subprocess.call(["rm", BAMunsorted])
    BAMsorted = os.path.abspath(BAMsorted)
    os.chdir("..")
    return BAMsorted


def _run_pileup(global_config, bamfile):
    """
    Perform samtools pileup on a .bam-file
    """

    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run samtools: samtools not present in "
                "the path and not in global config, please make sure to "
                "install samtools properly")

    pileupfile = bamfile.replace('.bam', '_coverage.csv')
    pileup_cmd = "{} mpileup {} | awk '{print $2, $4}' > {}".format(samtools,
            bamfile, pileupfile)
    p1 = subprocess.Popen(pileup_cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
    p1.wait()
    if p1.returncode == 0:
        return pileupfile
    else:
        print "Could not perform mpileup"
        return 1


def plot_coverage(pileupfile, samplename=None):
    """
    Plot coverage from samtools pileup output
    """
    if not samplename:
        samplename = pileupfile.split('.')[0]

    df = pd.io.parsers.read_csv(pileupfile, sep=' ', names=['pos', 'cov'])
    pl = df.plot(x='pos', y='cov')
    avg_cov = sum(df['cov'])/len(df)        #Calculate average coverage
    pl2 = plt.plot(range(len(df)), [avg_cov]*len(df), 'r--', linewidth=2)
    plt.xlim([0,len(df)])
    plt.ylabel('Coverage')
    plt.xlabel('Position')
    plt.title(samplename)
    plt.savefig(samplename)     #Save plot as .png
    return 0

