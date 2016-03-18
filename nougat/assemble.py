from __future__ import absolute_import
from __future__ import print_function
import sys, os, yaml, glob
import subprocess
import string
import sys
from nougat import common


def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    #Check if the user has specified tools, if not select default list of tools
    if "tools" not in sample_config or len(sample_config["tools"]) == 0:
        sample_config["tools"] = ["soapdenovo"]
    #Execute the commands now
    for command in sample_config["tools"]:
        command_fn = getattr( sys.modules[__name__] ,
                "_run_{}".format(command))
        sample_config = command_fn(global_config, sample_config,
                sorted_libraries_by_insert)


def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler                  = "abyss"
    outputName                 = sample_config["output"]
    currentDirectory           = os.getcwd()
    assemblyDirectory          = os.path.join(currentDirectory, assembler)
    # in abyss case there is no exectuable
    programBIN                 = global_config["Tools"][assembler]["bin"]
    program_options            = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure("abyss", assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART
    assembler_stdOut = open("abyss.stdOut", "a")
    assembler_stdErr = open("abyss.stdErr", "a")
    program=os.path.join(programBIN, "abyss-pe")

    command = ""
    command += "{} ".format(program)
    threads = 8 # default for UPPMAX
    if "threads" in sample_config :
        threads = sample_config["threads"]
    command += "np={} ".format(threads)
    kmer = 54
    if "kmer" in sample_config:
        kmer = sample_config["kmer"]
    command += "k={} ".format(kmer)

    libraries = {}
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                # check if this is the first time I insert a se file
                if "se" not in libraries:
                    libraries["se"] = "se=\'"
                libraries["se"] = libraries["se"] + read1
            else:
                if not "lib" in libraries:
                    libraries["lib"] = {}
                libName = insert # lib name is the insert size
                if not libName in libraries["lib"]:
                    libraries["lib"][libName] = ""
                libraries["lib"][libName] +=  "{} {} ".format(read1, read2)
        else:
            if not "mp" in libraries:
                libraries["mp"] = {}
            libName = format(insert)
            if not libName in libraries["mp"]:
                libraries["mp"][libName] = ""
            libraries["mp"][libName] +=  "{} {} ".format(read1, read2)
    #now create the command
    command += "name={} ".format(outputName)
    librariesSE = ""
    librariesPE = ""
    librariesMP = ""
    if "se" in libraries:
        libraries["se"] = libraries["se"] + "\'"
        librariesSE = libraries["se"]
    if "lib" in libraries:
        lib="lib=\'"
        for libPE, libPEreads in sorted(libraries["lib"].items()):
            lib = lib + "lib{} ".format(libPE)
            librariesPE += " lib{}=\'{}\' ".format(libPE,libPEreads)
        lib=lib + "\' "
        command += "{} ".format(lib)
    if "mp" in libraries:
        mp="mp=\'"
        for libMP, libMPreads in sorted(libraries["mp"].items()):
            mp = mp + "lib{} ".format(libMP)
            librariesMP += " lib{}=\'{}\' ".format(libMP,libMPreads)
        mp=mp + "\' "
        command += "{} ".format(mp)

    command += "{} ".format(librariesSE)
    command += "{} ".format(librariesPE)
    command += "{} ".format(librariesMP)

    common.print_command(command)
    if common.check_dryrun(sample_config):
        os.chdir("..")
        return sample_config

    os.makedirs(os.path.join(assemblyDirectory, "runABySS"))
    os.chdir("runABySS")
    returnValue = 0
    returnValue = subprocess.call(command, stdout=assembler_stdOut,
            stderr=assembler_stdErr, shell=True)
    os.chdir("..")
    flags = sample_config.get("flags", [])
    if returnValue == 0 and not common.check_dryrun(sample_config):
        if os.path.exists(os.path.join("runABySS","{}-contigs.fa".format(
            outputName))):
            subprocess.call(["cp", os.path.join("runABySS",
                "{}-contigs.fa".format(outputName)),
                "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join("runABySS",
                "{}-scaffolds.fa".format(outputName)),
                "{}.scf.fasta".format(outputName) ])
            if not "keep_tmp_files" in flags:
                subprocess.call(["rm", "-r", "runABySS"])
        elif not common.check_dryrun(sample_config):
            print("something wrong with ABySS -> no contig file generated")
            return sample_config
    else:
        print("ABySS terminated with an error. Please check running folder",
                "for more informations")
    os.chdir("..")
    return sample_config


def _run_allpaths(global_config, sample_config, sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler                  = "allpaths"
    outputName                 = sample_config["output"]
    currentDirectory           = os.getcwd()
    assemblyDirectory          = os.path.join(currentDirectory, assembler)
    # in abyss case there is no exectuable
    programBIN                 = global_config["Tools"][assembler]["bin"]
    program_options            = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure("allpaths", assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    inGroups_file = open("in_groups.csv", "w")
    inLibs_file   = open("in_libs.csv", "w")
    inGroups_file.write("group_name, library_name, file_name\n")
    inLibs_file.write("library_name, project_name, organism_name, type, "
            "paired, frag_size, frag_stddev, insert_size, insert_stddev, "
            "read_orientation,genomic_start, genomic_end\n")
    librariesForInLibs     = []
    librariesForInLibsDict = {}
    group_name             = 1;
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie":
            path, fqfile=os.path.split(read1)
            if "_1.fastq" in fqfile:
                fqfile = fqfile.replace("_1.fastq", "_?.fastq")
            elif "_R1_" in fqfile:
                fqfile = fqfile.replace("_R1_", "_R?_")
            else:
                print("error file format not supported {}".format(fqfile))
                return sample_config
            inGroups_file.write("PE{}, lib{}, {}\n".format(group_name, insert,
                os.path.join(path, fqfile)))
            group_name += 1
            if insert not in librariesForInLibsDict:
                librariesForInLibsDict[insert] = insert
                librariesForInLibs.append("lib{}, genome, genome, fragment, 1, "
                        "{}, {}, , , inward, 0, 0\n".format(insert,insert, std))
        elif orientation=="outtie":
            path, fqfile = os.path.split(read1)
            if "_1.fastq" in fqfile:
                fqfile = fqfile.replace("_1.fastq", "_?.fastq")
            elif "_R1_" in fqfile:
                fqfile = fqfile.replace("_R1_", "_R?_")
            else:
                print("error file format not supported {}".format(file))
                return sample_config
            inGroups_file.write("MP{}, lib{}, {}\n".format(group_name, insert,
                os.path.join(path, fqfile)))
            group_name += 1
            if insert not in librariesForInLibsDict:
                librariesForInLibsDict[insert] = insert
                librariesForInLibs.append("lib{}, genome, genome, fragment, 1, "
                        ", , {}, {}, outward, 0, 0\n".format(insert,insert, std))
        else:
            print("all paths support only innies and outties")
    inGroups_file.close()
    for lib in librariesForInLibs:
        inLibs_file.write(lib)
    inLibs_file.close()
    #NOW RUN ALLPATHS FOR REAL
    program=os.path.join(programBIN, "PrepareAllPathsInputs.pl")
    os.mkdir("data_dir")
    data_dir = os.path.join(os.getcwd(), "data_dir")
    ploidy = "PLOIDY=1"
    if len(program_options) > 0:
        if len(program_options) >1:
            print("Running ALlpaths only one parameter accepted as option",
                    "here: PLOIDY=2")
            return sample_config
        if program_options[0] == "PLOIDY=2":
            ploidy = "PLOIDY=2"
        else:
            print("Running ALlpaths only one parameter accepted as option",
                    "here: PLOIDY=2")
            return sample_config

    command = [program , "DATA_DIR={}".format(data_dir), ploidy, 
            "PICARD_TOOLS_DIR={}".format(
            global_config["Tools"]["picard"]["bin"]), 
            "FORCE_PHRED=True", "PHRED_64=False"]
    if common.check_dryrun(sample_config):
        common.print_command(command)
        program = os.path.join(programBIN, "RunAllPathsLG")
        command = [program, "PRE=.", "REFERENCE_NAME=.", "DATA_SUBDIR=data_dir",
                "RUN=allpaths", "SUBDIR=run"]
        common.print_command(command)
        os.chdir("..")
        return sample_config
    assembler_stdOut = open("allpaths_PrepareAllPathsInputs.stdOut", "w")
    assembler_stdErr = open("allpaths_PrepareAllPathsInputs.stdErr", "w")
    common.print_command(command)
    returnValue = subprocess.call(command,  stdout=assembler_stdOut, 
            stderr=assembler_stdErr)
    assembler_stdOut.close()
    assembler_stdErr.close()
    flags = sample_config.get("flags", [])
    if returnValue == 0:
        program = os.path.join(programBIN, "RunAllPathsLG")
        command = [program, "PRE=.", "REFERENCE_NAME=.", "DATA_SUBDIR=data_dir",
                "RUN=allpaths", "SUBDIR=run", "HAPLOIDIFY=True"]
        common.print_command(command)
        assembler_stdOut = open("allpaths_RunAllPathsLG.stdOut", "w")
        assembler_stdErr = open("allpaths_RunAllPathsLG.stdErr", "w")
        returnValue = subprocess.call(command,  stdout=assembler_stdOut,
                stderr=assembler_stdErr)
        if returnValue != 0:
            print("ALLPATHS RunAllPathsLG terminated with an error. Please",
                    "check running folder for more informations")
            os.chdir("..")
            return sample_config
        else: # save results
            assembly_dir = os.path.join("data_dir", "allpaths", "ASSEMBLIES",
                    "run")
            if os.path.exists(os.path.join(assembly_dir,
                "final.assembly.fasta")):
                exit_code = subprocess.call(["cp", os.path.join(assembly_dir,
                    "final.contigs.fasta"), "{}.ctg.fasta".format(outputName)])
                exit_code += subprocess.call(["cp", os.path.join(assembly_dir,
                    "final.assembly.fasta"), "{}.scf.fasta".format(outputName)])
                if not "keep_tmp_files" in flags and exit_code == 0:
                    subprocess.call(["rm", "-r", "data_dir"])
            else:
                print("something wrong with Allpaths > no contig file generated")
                os.chdir("..")
                return sample_config
    else:
        print("ALLPATHS PrepareAllPathInputs terminated with an error. "
                "Please check running folder for more informations")
        os.chdir("..")
        return sample_config
    os.chdir("..")
    return sample_config


def _run_cabog(global_config, sample_config, sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler = "cabog"
    outputName = sample_config["output"]
    currentDirectory = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    # in cabog case there is no exectuable
    programBIN = global_config["Tools"][assembler]["bin"]
    program_options = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure(assembler, assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART
    sys.path.insert(0, programBIN)
    libraries = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        command_fastqToCA = os.path.join(programBIN, "fastqToCA")
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        command_fastqToCA += " -libraryname "
        command_fastqToCA += " {}_{}".format(outputName, libraries)
        command_fastqToCA += " -insertsize "
        command_fastqToCA += " {} {} ".format(insert,std)
        command_fastqToCA += " -technology "
        command_fastqToCA += " illumina "
        command_fastqToCA += " -type "
        command_fastqToCA += " illumina "
        if orientation=="innie" or orientation=="none" :
            command_fastqToCA += " -innie "
            if read2 is None:
                command_fastqToCA += " -reads "
                command_fastqToCA += " {} ".format(read1)
            else:
                command_fastqToCA += " -mates "
                command_fastqToCA += " {},{} ".format(read1, read2)
        elif orientation=="outtie":
            command_fastqToCA += " -outtie "
            command_fastqToCA += " -mates "
            command_fastqToCA += " {},{} ".format(read1, read2)
        command_fastqToCA += " > "
        command_fastqToCA += " {}_{}.frg ".format(outputName, libraries)

        common.print_command(command_fastqToCA)
        if not common.check_dryrun(sample_config):
            cabog_stdOut = open("cabog_fastqToCA.stdOut", "w")
            cabog_stdErr = open("cabogfastqToCA.stdErr", "w")
            subprocess.call(command_fastqToCA, stderr=cabog_stdErr, shell=True)
            cabog_stdOut.close()
            cabog_stdErr.close()
        libraries += 1
    command_runCA = os.path.join(programBIN, "runCA")
    command_runCA += "  -d runCABOGfolder -p {} *frg".format(outputName)
    common.print_command(command_runCA)
    if common.check_dryrun(sample_config):
        return sample_config
    returnValue = 0
    cabog_stdOut = open("cabog_runCA.stdOut", "w")
    cabog_stdErr = open("cabog_runCA.stdErr", "w")
    returnValue = subprocess.call(command_runCA, stdout=cabog_stdOut,
            stderr=cabog_stdErr, shell=True)
    flags = sample_config.get("flags", [])
    if returnValue == 0:
        #assembly succed, remove files and save assembly
        if os.path.exists(os.path.join("runCABOGfolder","9-terminator",
            "{}.ctg.fasta".format(outputName))):
            subprocess.call(["cp", os.path.join("runCABOGfolder","9-terminator",
                "{}.ctg.fasta".format(outputName)),
                "{}.ctg.fasta".format(outputName)])
            subprocess.call(["cp", os.path.join("runCABOGfolder","9-terminator",
                "{}.scf.fasta".format(outputName)),
                "{}.scf.fasta".format(outputName)])
            if not "keep_tmp_files" in flags:
                subprocess.call(["rm", "-r", "runCABOGfolder"])
        else:
            print("something wrong with CABOG -> no contig file generated")
    else:
        print("CABOG terminated with an error. Please check running folder",
                "for more informations")
    os.chdir("..")
    return sample_config


def _run_masurca(global_config, sample_config,sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler = "masurca"
    outputName = sample_config["output"]
    currentDirectory = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    # in cabog case there is no exectuable
    programBIN = global_config["Tools"][assembler]["bin"]
    program_options = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure(assembler, assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    masurca_config_file = open("configuration.txt", "w")
    masurca_config_file.write("DATA\n")
    allTheLetters = string.lowercase
    libraryPE    = "p"
    libraryPEnum = 0
    libraryMP    = "m"
    libraryMPnum = 0
    #TODO: single ended reads
    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        if orientation=="innie":
            if read2 is not None:
                configurationLine = "PE = {}{} {} {} {} {}".format(libraryPE,
                        allTheLetters[libraryPEnum], insert, std, read1, read2)
                masurca_config_file.write("{}\n".format(configurationLine))
                libraryPEnum += 1
                #TODO: check when more than 21 PE libraries ae specified
        elif orientation=="outtie":
            configurationLine = "JUMP = {}{} {} {} {} {}".format(libraryMP, 
                    allTheLetters[libraryMPnum], insert, std, read1, read2)
            masurca_config_file.write("{}\n".format(configurationLine))
            libraryMPnum += 1
            #TODO: check when more than 21 PE libraries ae specified
    masurca_config_file.write("END\n")

    masurca_config_file.write("\n")

    masurca_config_file.write("PARAMETERS\n")
    #this is k-mer size for deBruijn graph values between 25 and 101 are 
    #supported, auto will compute the optimal size based on the read data 
    #and GC content
    masurca_config_file.write("GRAPH_KMER_SIZE=auto\n")
    #set this to 1 for Illumina-only assemblies and to 0 if you have 2x or 
    #more long (Sanger, 454) reads
    masurca_config_file.write("USE_LINKING_MATES=1\n")
    #this parameter is useful if you have too many jumping library mates. 
    #See manual for explanation about settings based on genome length
    if sample_config["genomeSize"] > 10000000:
        masurca_config_file.write("LIMIT_JUMP_COVERAGE = 1000\n")
    else:
        masurca_config_file.write("LIMIT_JUMP_COVERAGE = 60\n")
    #these are the additional parameters to Celera Assembler.  do not worry 
    #about performance, number or processors or batch sizes -- these are 
    #computed automatically. for mammals do not set cgwErrorRate above 0.15!!!
    if sample_config["genomeSize"] > 1500000000:
        masurca_config_file.write("CA_PARAMETERS = ovlMerSize=30 \
                cgwErrorRate=0.15 ovlMemory=4GB\n")
    else:
        masurca_config_file.write("CA_PARAMETERS = ovlMerSize=30 \
                cgwErrorRate=0.25 ovlMemory=4GB\n")
    #auto-detected number of cpus to use
    threads = 8 # default for UPPMAX
    if "threads" in sample_config :
        threads = sample_config["threads"]
    masurca_config_file.write("NUM_THREADS= {}\n".format(threads))
    #this is mandatory jellyfish hash size ---- jellyfish hash size, 
    #set this to about 10x the genome size.
    JF_SIZE = sample_config["genomeSize"] * 11
    masurca_config_file.write("JF_SIZE={}\n".format(JF_SIZE))
    #this specifies if we do (1) or do not (0) want to trim long runs of 
    #homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
    masurca_config_file.write("DO_HOMOPOLYMER_TRIM=0\n")
    masurca_config_file.write("END\n")
    masurca_config_file.write("\n")

    masurca_config_file.close()

    if common.check_dryrun(sample_config):
        os.chdir("..")
        return sample_config

    masurca_stdOut = open("masurca.stdOut", "w")
    masurca_stdErr = open("masurca.stdErr", "w")
    os.mkdir("runMASURCA")
    os.chdir("runMASURCA")
    command = [os.path.join(programBIN,"bin/masurca") , "../configuration.txt"]
    common.print_command(command)

    subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)
    if not os.path.exists("assemble.sh"):
        print("MaSuRCA: assemble.sh not created. Unknown failure")
        return sample_config
    command = ["./assemble.sh"]
    common.print_command(command)
    returnValue = subprocess.call(command, stdout=masurca_stdOut,
            stderr=masurca_stdErr)
    os.chdir("..")
    flags = sample_config.get("flags", [])
    if returnValue == 0:
        if os.path.exists(os.path.join(
            "runMASURCA","CA/10-gapclose/genome.scf.fasta")):
            subprocess.call(["cp", os.path.join(
                "runMASURCA","CA/10-gapclose/genome.ctg.fasta"), 
                "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join(
                "runMASURCA","CA/10-gapclose/genome.scf.fasta"),
                "{}.scf.fasta".format(outputName) ])
            if not "keep_tmp_files" in flags:
                subprocess.call(["rm", "-r", "runMASURCA"])
        else:
            print("something wrong with MaSuRCA -> no contig file generated")
    else:
        print("MaSuRCA terminated with an error. Please check running folder",
                "for more informations")
        return sample_config
    os.chdir("..")
    return sample_config


def _run_soapdenovo(global_config, sample_config, sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler = "soapdenovo"
    outputName = sample_config["output"]
    currentDirectory = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    # in cabog case there is no exectuable
    programBIN = global_config["Tools"][assembler]["bin"]
    program_options = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure(assembler, assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART
    kmer = 54
    if "kmer" in sample_config:
        kmer = sample_config["kmer"]
    threads = ["-p", "8"] # default for UPPMAX
    if "threads" in sample_config:
        threads = ["-p", "{}".format(sample_config["threads"])]
    soap_config_file = open("configuration.txt", "w")
    soap_config_file.write("max_rd_len=100\n")
    #TODO make this a parameter in the options
    rank = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        soap_config_file.write("[LIB]\n")
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        soap_config_file.write("avg_ins={}\n".format(insert))
        soap_config_file.write("rank={}\n".format(rank))
        rank += 1
        soap_config_file.write("map_len=30\n")
        if orientation=="innie" or orientation=="none":
            soap_config_file.write("asm_flags=3\n")
            soap_config_file.write("pair_num_cutoff=3\n")
            soap_config_file.write("reverse_seq=0\n")
            if read2 is None:
                soap_config_file.write("q={}\n".format(read1))
            else:
                soap_config_file.write("q1={}\n".format(read1))
                soap_config_file.write("q2={}\n".format(read2))
        elif orientation=="outtie":
            soap_config_file.write("asm_flags=2\n")
            soap_config_file.write("pair_num_cutoff=5\n")
            soap_config_file.write("reverse_seq=1\n")
            soap_config_file.write("q1={}\n".format(read1))
            soap_config_file.write("q2={}\n".format(read2))

    soap_config_file.close()
    assembler_stdOut = open("soap.stdOut", "w")
    assembler_stdErr = open("soap.stdErr", "w")
    os.makedirs(os.path.join(assemblyDirectory, "runSOAP"))
    os.chdir("runSOAP")
    #TODO : lots of missing options
    command = [programBIN , "all", "-s", "../configuration.txt", "-K",
            "{}".format(kmer), "-L", "500", "-o", "soapAssembly", threads[0],
            threads[1] ]
    common.print_command(command)
    returnValue = 0
    if not common.check_dryrun(sample_config):
        subprocess.call(command, stdout=assembler_stdOut,
                stderr=assembler_stdErr)
    else:
        os.chdir("..")
        os.chdir("..")
        return sample_config

    os.chdir("..")
    flags = sample_config.get("flags", [])
    if returnValue == 0:
        if(os.path.exists(os.path.join("runSOAP","soapAssembly.scafSeq"))):
            subprocess.call(["cp", os.path.join("runSOAP",
                "soapAssembly.scafSeq"), "{}.scf.fasta".format(outputName)])
            subprocess.call(["cp", os.path.join("runSOAP",
                "soapAssembly.contig"), "{}.ctg.fasta".format(outputName)])
            if not "keep_tmp_files" in flags:
                subprocess.call(["rm", "-r", "runSOAP"])
        else:
            print("something wrong with SOAPdenovo -> no contig file generated")
    else:
        print("SOAPdenovo terminated with an error. Please check running",
                "folder for more informations")
        os.chdir("..")
        return sample_config
    os.chdir("..")
    return sample_config


def _run_spades(global_config, sample_config, sorted_libraries_by_insert):
    ########## ACQUIRE ALL THE INFO AND CREATE THE ASSEMBLY FOLDER
    assembler = "spades"
    outputName = sample_config["output"]
    currentDirectory = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    # in cabog case there is no exectuable
    programBIN = global_config["Tools"][assembler]["bin"]
    program_options = global_config["Tools"][assembler]["options"]
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if _prepare_folder_structure(assembler, assemblyDirectory) == 0:
        os.chdir(assemblyDirectory)
    else:
        return sample_config
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    command = ""
    command += "{} ".format(programBIN)
    for option in program_options:
        command += "{} ".format(option)

    #creates the command on-the-fly
    peLibrary = 1
    mpLibrary = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                command += "--pe{}-s {} ".format(peLibrary, read1)
            else:
                command += "--pe{}-1 {} --pe{}-2 {} ".format(peLibrary, read1,
                        peLibrary, read2)
            peLibrary += 1
        elif orientation=="outtie":
            command += "--mp{}-1 {} --mp{}-2 {} ".format(mpLibrary, read1,
                    mpLibrary, read2)
            mpLibrary += 1
        else:
            print("orientation{} not supported.... why the program did not",
                    "failed earlier?".format(orientation))

    command += "-o {} ".format(outputName)
    common.print_command(command)
    returnValue = 0
    if not common.check_dryrun(sample_config):
        assembler_stdOut = open("spades.stdOut", "a")
        assembler_stdErr = open("spades.stdErr", "a")
        returnValue = subprocess.call(command, stdout=assembler_stdOut,
                stderr=assembler_stdErr, shell=True)
    else:
        return sample_config

    flags = sample_config.get("flags", [])
    if returnValue == 0:
        if os.path.exists(os.path.join(outputName,"contigs.fasta")):
            subprocess.call(["cp", os.path.join(outputName,"contigs.fasta"),
                "{}.ctg.fasta".format(outputName)])
            subprocess.call(["cp", os.path.join(outputName,"scaffolds.fasta"),
                "{}.scf.fasta".format(outputName)])
            if not "keep_tmp_files" in flags:
                subprocess.call(["rm", "-r", outputName])
        else:
            print("something wrong with SPADES -> no contig file generated")
    else:
        print("SPADES terminated with an error. Please check running folder",
                "for more informations")

    os.chdir("..")
    return sample_config


def _run_trinity(global_config, sample_config, sorted_libraries_by_insert):
    print("running trinity ...")
    assembler = "trinity"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(
            sorted_libraries_by_insert)
     # in masurca case there is no exectuable as a make file must be created
    programBIN = global_config["Tools"][assembler]["bin"] + "Trinity"
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    command = [programBIN]
    command.extend(["--seqType", "fq"])
    command.extend(["--JM", "100G"])
    if "threads" in sample_config:
        command.extend(["--CPU", str(sample_config["threads"])])

    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        if read2 is None:
            command.append("--single")
            command.append("{}".format(read1))
        elif orientation=="innie":
            command.append("--left")
            command.append("{}".format(read1))
            command.append("--right")
            command.append("{}".format(read2))
        else:
            print("trinity: somthing wrong or unexpected in the sample",
                    "config file")
            return sample_config
    command.extend(["--output", "trinity"])
    assembler_stdOut = open("trinity.stdOut", "w")
    assembler_stdErr = open("trinity.stdErr", "w")
    print(" ".join(command))

    returnValue = subprocess.call(" ".join(command), stdout=assembler_stdOut,
            stderr=assembler_stdErr, shell=True)

    # now align reads back to transcripts and estimate abundance
    os.chdir("trinity")
    programBIN = global_config["Tools"][assembler]["bin"] + \
            "util/align_and_estimate_abundance.pl"
    command = [programBIN]
    command.extend(["--transcripts", "Trinity.fasta"])
    command.extend(["--seqType", "fq"])
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        if read2 is not None and orientation == "innie":
            command.append("--left")
            command.append("{}".format(os.path.splitext(read1)[0]))
            command.append("--right")
            command.append("{}".format(os.path.splitext(read2)[0]))

    command.extend(["--aln_method", "bowtie"])
    command.extend(["--est_method", "RSEM"])
    command.append("--debug")
    command.append("--trinity_mode")
    command.append("--prep_reference")

    if "threads" in sample_config:
        command.extend(["--thread_count", str(sample_config["threads"])])
    print(" ".join(command))
    returnValue = subprocess.call(" ".join(command), stdout=assembler_stdOut,
            stderr=assembler_stdErr, shell=True)

    #now copy results
    os.chdir("..")
    subprocess.call(["cp", "trinity/Trinity.fasta",
        "{}.fasta".format(outputName)])
    subprocess.call(["cp", "trinity/RSEM.isoforms.results",
        "{}.isoforms.results".format(outputName)])
    subprocess.call(["cp", "trinity/RSEM.genes.results",
        "{}.genes.results".format(outputName)])
    os.chdir(currentDirectory)
    return sample_config


def _prepare_folder_structure(assembler,assemblyDirectory):
    if common.directory_exists(assemblyDirectory):
        print("Assembler {} asumer already computed as folder {} exists".format(
            assembler,assemblyDirectory))
        return 1
    return 0


def _run_abyss_mergePairs(global_config, sample_config, 
        sorted_libraries_by_insert):
    print("running abyss-mergepairs ...")
    assembler = "abyss_mergePairs"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(
            sorted_libraries_by_insert)
    # in abyss case there is no exectuable
    programBIN      = global_config["Tools"][assembler]["bin"]
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    program=programBIN
    command = []
    command.append(program)
    for option in program_options:
            command.append(option)

    libraries = {}
    for library, libraryInfo in sorted_libraries_by_insert:
        read1 = libraryInfo["pair1"]
        read2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        std = libraryInfo["std"]
        outputNameArray  = read1.split('/')[-1].split('_')
        outputName = "{}_{}".format(outputNameArray[0], outputNameArray[1])

        if orientation=="innie":
            if read2 is not None:
                currentCommand = command;
                currentCommand.append('-o')
                currentCommand.append(outputName)
                currentCommand.append(read1)
                currentCommand.append(read2)
                abyss_stdOut = open("mergePairs_{}.stdOut".format(outputName),
                        "a")
                abyss_stdErr = open("mergePairs_{}.stdErr".format(outputName),
                        "a")
                print(command)
                subprocess.call(command, stdout=abyss_stdOut,
                        stderr=abyss_stdErr)
                command_mv = ["mv", "mergePairs_{}.stdErr".format(outputName),
                        "{}.txt".format(outputName)]
                subprocess.call(command_mv)

    os.chdir("..")
    return sample_config

