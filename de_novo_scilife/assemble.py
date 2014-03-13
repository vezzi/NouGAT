import sys, os, yaml, glob
import subprocess
import string
import sys
import common




def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if "tools" in sample_config:
        #need to follow the commands listed here
        for command in sample_config["tools"]:
            command_fn = getattr( sys.modules[__name__] , "_run_{}".format(command))
            sample_config = command_fn(global_config, sample_config, sorted_libraries_by_insert)
    else:
        #run default pipeline for de novo assembly
        sample_config = _run_soapdenovo(global_config, sample_config, sorted_libraries_by_insert)
    

def _run_trinity(global_config, sample_config, sorted_libraries_by_insert):
    print "running trinity ..."
    assembler = "trinity"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] + "Trinity.pl"  # in masurca case there is no exectuable as a make file must be created
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    command = [programBIN]
    command.append("--seqType")
    command.append("fq")
    command.append("--JM")
    command.append("50G")
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if read2 is None:
            command.append("--single")
            command.append("{}".format(read1))
        elif orientation=="innie":
            command.append("--left")
            command.append("{}".format(read1))
            command.append("--right")
            command.append("{}".format(read2))
        else:
            print "trinity: somthing wrong or unexpected in the sample config file"
            return sample_config
    command.append("--output")
    command.append("trinity")
    assembler_stdOut = open("trinity.stdOut", "w")
    assembler_stdErr = open("trinity.stdErr", "w")
    print command

    returnValue = subprocess.call(" ".join(command), stdout=assembler_stdOut, stderr=assembler_stdErr, shell=True)

    # now align reads back to transcripts
    os.chdir("trinity")
    programBIN = global_config["Tools"][assembler]["bin"] + "util/alignReads.pl"
    command = [programBIN]
    command.append("--target")
    command.append("Trinity.fasta")
    command.append("--seqType")
    command.append("fq")
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if read2 is not None and orientation == "innie":
            command.append("--left")
            command.append("{}".format(os.path.splitext(read1)[0]))
            command.append("--right")
            command.append("{}".format(os.path.splitext(read2)[0]))

    command.append("--aligner")
    command.append("bowtie")
    command.append("--retain_intermediate_files")
    print command
    returnValue = subprocess.call(" ".join(command), stdout=assembler_stdOut, stderr=assembler_stdErr, shell=True)

    # now quantify trnascripts
    programBIN = global_config["Tools"][assembler]["bin"] + "util/RSEM_util/run_RSEM_align_n_estimate.pl"
    command = [programBIN]
    command.append("--transcripts")
    command.append("Trinity.fasta")
    command.append("--seqType")
    command.append("fq")
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if read2 is not None and orientation == "innie":
            command.append("--left")
            command.append("{}".format(os.path.splitext(read1)[0]))
            command.append("--right")
            command.append("{}".format(os.path.splitext(read2)[0]))

    print command
    returnValue = subprocess.call(" ".join(command), stdout=assembler_stdOut, stderr=assembler_stdErr, shell=True)

    #now copy results
    os.chdir("..")
    subprocess.call(["cp", "trinity/Trinity.fasta", "{}.fasta".format(outputName)])
    subprocess.call(["cp", "trinity/RSEM.isoforms.results", "{}.isoforms.results".format(outputName)])
    subprocess.call(["cp", "trinity/RSEM.genes.results", "{}.genes.results".format(outputName)])
    os.chdir(currentDirectory)
    return sample_config




def _run_allpaths(global_config, sample_config, sorted_libraries_by_insert):
    print "running allpaths ..."
    assembler = "allpaths"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART

    inGroups_file = open("in_groups.csv", "w")
    inLibs_file   = open("in_libs.csv", "w")
    inGroups_file.write("group_name, library_name, file_name\n")
    inLibs_file.write("library_name, project_name, organism_name, type, paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation,genomic_start, genomic_end\n")
    
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
            path,file=os.path.split(read1)
            file = file.replace("_1.fastq", "_?.fastq")
            inGroups_file.write("PE{}, lib{}, {}\n".format(group_name, insert, os.path.join(path, file)))
            group_name += 1
            if insert not in librariesForInLibsDict:
                librariesForInLibsDict[insert] = insert
                librariesForInLibs.append("lib{}, genome, genome, fragment, 1, {}, {}, , , inward, 0, 0\n".format(insert,insert, std))
        elif orientation=="outtie":
            path,file=os.path.split(read1)
            file = file.replace("_1.fastq", "_?.fastq")
            inGroups_file.write("MP{}, lib{}, {}\n".format(group_name, insert, os.path.join(path, file)))
            group_name += 1
            if insert not in librariesForInLibsDict:
                librariesForInLibsDict[insert] = insert
                librariesForInLibs.append("lib{}, genome, genome, fragment, 1, , , {}, {}, outward, 0, 0\n".format(insert,insert, std))
        else:
            print "all paths support only innies and outties"

    inGroups_file.close()
    for lib in librariesForInLibs:
        inLibs_file.write(lib)
    inLibs_file.close()

    #NOW RUN ALLPATHS FOR REAL
    #TODO: must check if path is correctly set
    ### PrepareAllPathsInputs.pl DATA_DIR=$PWD PLOIDY=1 PICARD_TOOLS_DIR=/home/francesco.vezzi/DE_NOVO_PIPELINE/tools/picard-tools-1.103/ FORCE_PHRED=True PHRED_64=False
    command = ["PrepareAllPathsInputs.pl" , "DATA_DIR=$PWD", "PLOIDY=1", "PICARD_TOOLS_DIR={}".format(global_config["Tools"]["picard"]["bin"]), "FORCE_PHRED=True", "PHRED_64=False"]
    print command
    returnValue = subprocess.call(command)
    if returnValue == 0:
        command = ["RunAllPathsLG", "PRE=.", "REFERENCE_NAME=.", "DATA_SUBDIR=.", "RUN=allpaths", "SUBDIR=run"]
        returnValue = subprocess.call(command)
        if returnValue != 0:
            print "ALLPATHS RunAllPathsLG terminated with an error. Please check running folder for more informations"
            oc.chdir("..")
            return sample_config
    else:
        print "ALLPATHS PrepareAllPathInputs terminated with an error. Please check running folder for more informations"
        oc.chdir("..")
        return sample_config
    os.chdir("..")
    return sample_config



def _run_soapdenovo(global_config, sample_config, sorted_libraries_by_insert):
    print "running SOAPdenovo ..."
    assembler = "soapdenovo"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    
    
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART
    kmer = 54
    if "kmer" in sample_config:
        kmer = sample_config["kmer"]
    threads = ["-p", "8"] # default for UPPMAX
    if "threads" in sample_config:
        threads = ["-p", "{}".format(sample_config["threads"])]

    soap_config_file = open("configuration.txt", "w")
    soap_config_file.write("max_rd_len=100\n") #TODO make this a parameter in the options
    rank = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        soap_config_file.write("[LIB]\n")
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
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
    command = [programBIN , "all", "-s", "../configuration.txt", "-K", "{}".format(kmer), "-L", "500", "-o", "soapAssembly", threads[0] , threads[1] ]
    print command
    


    returnValue = subprocess.call(command, stdout=assembler_stdOut, stderr=assembler_stdErr)

    os.chdir("..")
    if returnValue == 0:
        if(os.path.exists(os.path.join("runSOAP","soapAssembly.scafSeq"))):
            subprocess.call(["mv", os.path.join("runSOAP","soapAssembly.scafSeq"), "{}.scf.fasta".format(outputName) ])
            subprocess.call(["mv", os.path.join("runSOAP","soapAssembly.contig"), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runSOAP"])
        else:
            print "something wrong with SOAPdenovo -> no contig file generated"
    else:
        print "SOAPdenovo terminated with an error. Please check running folder for more informations"
        oc.chdir("..")
        return sample_config
    os.chdir("..")
    return sample_config


def _run_masurca(global_config, sample_config,sorted_libraries_by_insert):
    print "running MaSuRCA ..."
    print "MaSurCA still to be implemented "
    return sample_config
    
    assembler = "masurca"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    
    
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART
    programPATH = global_config["assemble"]["masurca"]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options=global_config["assemble"]["masurca"]["options"]

    masurca_config_file = open("configuration.txt", "w")
    masurca_config_file.write("PATHS\n")
    masurca_config_file.write("JELLYFISH_PATH=" + os.path.join("{}".format(programPATH), "bin" ) + "\n" )
    masurca_config_file.write("SR_PATH=" +  os.path.join("{}".format(programPATH), "bin" ) + "\n" )
    masurca_config_file.write("CA_PATH=" +  os.path.join("{}".format(programPATH), "CA/Linux-amd64/bin") + "\n")
    masurca_config_file.write("END\n")
    
    masurca_config_file.write("\n")
    
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
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie":
            if read2 is not None:
                configurationLine = "PE = {}{} {} {} {} {}".format(libraryPE, allTheLetters[libraryPEnum], insert, std, read1, read2)
                masurca_config_file.write("{}\n".format(configurationLine))
                libraryPEnum+=1 #TODO: check when more than 21 PE libraries ae specified
        elif orientation=="outtie":
            configurationLine = "JUMP = {}{} {} {} {} {}".format(libraryMP, allTheLetters[libraryMPnum], insert, std, read1, read2)
            masurca_config_file.write("{}\n".format(configurationLine))
            libraryMPnum += 1  #TODO: check when more than 21 PE libraries ae specified
    masurca_config_file.write("END\n")
    
    masurca_config_file.write("\n")
    
    masurca_config_file.write("PARAMETERS\n")
    #this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content
    masurca_config_file.write("GRAPH_KMER_SIZE=auto\n")
    #set this to 1 for Illumina-only assemblies and to 0 if you have 2x or more long (Sanger, 454) reads
    masurca_config_file.write("USE_LINKING_MATES=1\n")
    #this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and something large (300) for mammals
    masurca_config_file.write("LIMIT_JUMP_COVERAGE = 60\n")
    #these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. for mammals do not set cgwErrorRate above 0.15!!!
    masurca_config_file.write("CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.25 ovlMemory=4GB\n")
    #minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100
    masurca_config_file.write("KMER_COUNT_THRESHOLD = 1\n")
    #auto-detected number of cpus to use
    masurca_config_file.write("NUM_THREADS= 8\n")
    #this is mandatory jellyfish hash size
    masurca_config_file.write("JF_SIZE=1000000000\n")
    #this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
    masurca_config_file.write("DO_HOMOPOLYMER_TRIM=0\n")
    masurca_config_file.write("END\n")
    masurca_config_file.write("\n")

    masurca_config_file.close()

    masurca_stdOut = open("masurca.stdOut", "w")
    masurca_stdErr = open("masurca.stdErr", "w")
    command = [os.path.join(programPATH,"bin/runSRCA.pl") , "configuration.txt"]
    subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)

    if not os.path.exists("assemble.sh"):
        print "MaSuRCA: assemble.sh not created. Unknown failure"
        return sample_config
    command = ["./assemble.sh"]
    returnValue = subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)
    if returnValue == 0:
        if os.path.exists(os.path.join("runABySS","{}-contigs.fa".format(outputName))):
            subprocess.call(["cp", os.path.join("runABySS","{}-contigs.fa".format(outputName)), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join("runABySS","{}-scaffolds.fa".format(outputName)), "{}.scaf.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runABySS"])
        else:
            print "something wrong with MaSuRCA -> no contig file generated"
    else:
        print "MaSuRCA terminated with an error. Please check running folder for more informations"
        return sample_config
    os.chdir("..")
    return sample_config





def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
    print "running ABySS ..."
    assembler = "abyss"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in abyss case there is no exectuable
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]

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
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                if "se" not in libraries: # check if this is the first time I insert a se file
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
    librariesSE       = ""
    librariesPE       = ""
    librariesMP       = ""
    if "se" in libraries:
        libraries["se"] = libraries["se"] + "\'"
        librariesSE = libraries["se"]
    if "lib" in libraries:
        lib="lib=\'"
        for libPE, libPEreads in sorted(libraries["lib"].iteritems()):
            lib = lib + "lib{} ".format(libPE)
            librariesPE += " lib{}=\'{}\' ".format(libPE,libPEreads)
        lib=lib + "\' "
        command += "{} ".format(lib)
    if "mp" in libraries:
        mp="mp=\'"
        for libMP, libMPreads in sorted(libraries["mp"].iteritems()):
            mp = mp + "lib{} ".format(libMP)
            librariesMP += " lib{}=\'{}\' ".format(libMP,libMPreads)
        mp=mp + "\' "
        command += "{} ".format(mp)

    command += "{} ".format(librariesSE)
    command += "{} ".format(librariesPE)
    command += "{} ".format(librariesMP)

    os.makedirs(os.path.join(assemblyDirectory, "runABySS"))
    os.chdir("runABySS")
    print command
    returnValue = subprocess.call(command, stdout=assembler_stdOut, stderr=assembler_stdErr, shell=True)

    os.chdir("..")
    if returnValue == 0:
        if os.path.exists(os.path.join("runABySS","{}-contigs.fa".format(outputName))):
            subprocess.call(["cp", os.path.join("runABySS","{}-contigs.fa".format(outputName)), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join("runABySS","{}-scaffolds.fa".format(outputName)), "{}.scf.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runABySS"])
        else:
            print "something wrong with ABySS -> no contig file generated"
            return sample_config
    else:
        print "ABySS terminated with an error. Please check running folder for more informations"
    os.chdir("..")
    return sample_config





def _run_abyss_mergePairs(global_config, sample_config , sorted_libraries_by_insert):
    print "running abyss-mergepairs ..."
    assembler = "abyss_mergePairs"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in abyss case there is no exectuable
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
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        outputNameArray  = read1.split('/')[-1].split('_')
        outputName = "{}_{}".format(outputNameArray[0], outputNameArray[1])
        
        if orientation=="innie":
            if read2 is not None:
                currentCommand = command;
                currentCommand.append('-o')
                currentCommand.append(outputName)
                currentCommand.append(read1)
                currentCommand.append(read2)
                abyss_stdOut = open("mergePairs_{}.stdOut".format(outputName), "a")
                abyss_stdErr = open("mergePairs_{}.stdErr".format(outputName), "a")
                print command
                subprocess.call(command, stdout=abyss_stdOut, stderr=abyss_stdErr)
                command_mv = ["mv", "mergePairs_{}.stdErr".format(outputName), "{}.txt".format(outputName)]
                subprocess.call(command_mv)

    os.chdir("..")
    return sample_config






def _run_spades(global_config, sample_config, sorted_libraries_by_insert):
    print "running spades ..."
    assembler = "spades"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in abyss case there is no exectuable
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART


    assembler_stdOut = open("spades.stdOut", "a")
    assembler_stdErr = open("spades.stdErr", "a")
    
    command = ""
    command += "{} ".format(programBIN)
    for option in program_options:
        command += "{} ".format(option)

    #creates the command on-the-fly
    peLibrary = 1
    mpLibrary = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                command += "--pe{}-s {} ".format(peLibrary, read1)
            else:
                command += "--pe{}-1 {} --pe{}-2 {} ".format(peLibrary, read1, peLibrary, read2)
            peLibrary += 1
        elif orientation=="outtie":
            command += "--mp{}-1 {} --mp{}-2 {} ".format(mpLibrary, read1, mpLibrary, read2)
            mpLibrary += 1
        else:
            print "orientation{} not supported.... why the program didnot failed earlier?".format(orientation)

    command += "-o {} ".format(outputName)

    returnValue = subprocess.call(command, stdout=assembler_stdOut, stderr=assembler_stdErr, shell=True)
    if returnValue == 0:
        if os.path.exists(os.path.join(outputName,"contigs.fasta")):
            subprocess.call(["cp", os.path.join(outputName,"contigs.fasta"),  "{}.ctg.fasta".format(outputName)])
            subprocess.call(["cp", os.path.join(outputName,"scaffolds.fasta"),  "{}.scf.fasta".format(outputName)])
            subprocess.call(["rm", "-r", outputName])
        else:
            print "something wrong with SPADES -> no contig file generated"
    else:
        print "SPADES terminated with an error. Please check running folder for more informations"

    os.chdir("..")
    return sample_config



def _run_cabog(global_config, sample_config, sorted_libraries_by_insert):
    print "running cabog ..."
    assembler = "cabog"
    outputName = sample_config["output"]
    currentDirectory  = os.getcwd()
    assemblyDirectory = os.path.join(currentDirectory, assembler)
    if common.directory_exists(assemblyDirectory):
        return sample_config
    os.chdir(assemblyDirectory) # now I am in the assembly directory
    sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)
    programBIN      = global_config["Tools"][assembler]["bin"] # in abyss case there is no exectuable
    program_options = global_config["Tools"][assembler]["options"]
    if assembler in sample_config:
        program_options=sample_config[assembler]
    ########### HERE IT START THE SPECIFIC ASSEMBLER PART


    sys.path.insert(0, programBIN)
    libraries = 1
    cabog_stdOut = open("cabog.stdOut", "w")
    cabog_stdErr = open("cabog.stdErr", "w")
    for library, libraryInfo in sorted_libraries_by_insert:
        command_fastqToCA = os.path.join(programBIN, "fastqToCA")
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        command_fastqToCA  += " -libraryname "
        command_fastqToCA  += " {}_{}".format(outputName, libraries)
        command_fastqToCA  += " -insertsize "
        command_fastqToCA  += " {} {} ".format(insert,std)
        command_fastqToCA  += " -technology "
        command_fastqToCA  += " illumina "
        command_fastqToCA  += " -type "
        command_fastqToCA  += " illumina "
        
        if orientation=="innie" or orientation=="none" :
            command_fastqToCA  += " -innie "
            if read2 is None:
                command_fastqToCA  += " -reads "
                command_fastqToCA  += " {} ".format(read1)
            else:
                command_fastqToCA  += " -mates "
                command_fastqToCA  += " {},{} ".format(read1, read2)
        elif orientation=="outtie":
            command_fastqToCA  += " -outtie "
            command_fastqToCA  += " -mates "
            command_fastqToCA  += " {},{} ".format(read1, read2)
        command_fastqToCA  += " > "
        command_fastqToCA  += " {}_{}.frg ".format(outputName, libraries)
        subprocess.call(command_fastqToCA, stderr=cabog_stdErr, shell=True)
        libraries += 1
    
    command_runCA = os.path.join(programBIN, "runCA")
    command_runCA += "  -d runCABOGfolder -p {} *frg".format(outputName)
    returnValue = subprocess.call(command_runCA, stdout=cabog_stdOut, stderr=cabog_stdErr, shell=True)
    if returnValue == 0:
        #assembly succed, remove files and save assembly
        if os.path.exists(os.path.join("runCABOGfolder","9-terminator", "{}.ctg.fasta".format(outputName))):
            subprocess.call(["mv", os.path.join("runCABOGfolder","9-terminator", "{}.ctg.fasta".format(outputName)), "{}.ctg.fasta".format(outputName)])
            subprocess.call(["mv", os.path.join("runCABOGfolder","9-terminator", "{}.scf.fasta".format(outputName)), "{}.scf.fasta".format(outputName)])
            subprocess.call(["rm", "-r", "runCABOGfolder"])
        else:
            print "something wrong with CABOG -> no contig file generated"
    else:
        print "CABOG terminated with an error. Please check running folder for more informations"

    os.chdir("..")
    return sample_config
