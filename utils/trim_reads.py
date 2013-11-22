import sys, os, yaml, glob
import subprocess
import argparse
from copy import copy, deepcopy



def main(args):
    workingDir = os.getcwd()
    samples_data_dir = args.sample_data_dir
        
    if(args.global_config is None):
        sys.exit("Error: global configuration must be speficied")
            
    with open(args.global_config) as in_handle:
        global_config = yaml.load(in_handle)
        
    command = ["java", "-jar"];
    if "trimmomatic" not in global_config["Tools"]:
        sys.exit("Error: trimmomatic is not specified in the global configuratipon yaml file")
    program=global_config["Tools"]["trimmomatic"]["bin"]
    program_options=global_config["Tools"]["trimmomatic"]["options"]

    command.append(program)
    for option in program_options:
        command.append(option)

    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        currentCommand = []
        currentCommand = deepcopy(command)
        SampleQC_folder = os.path.join(os.getcwd(), sample_dir_name)
        print sample_dir_name
        if not os.path.exists(SampleQC_folder):
            os.makedirs(SampleQC_folder)
        else:
            print "done ({} folder already present, assumed already run)".format(SampleQC_folder)
            continue

        os.chdir(SampleQC_folder)
        sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
        sample_files = [ f for f in os.listdir(sample_data_dir) if os.path.isfile(os.path.join(sample_data_dir,f))]
        
        pair1_file = "";
        pair2_file = "";
        for file in sample_files:
            if "_R1_" in file:
                pair1_file = os.path.join(sample_data_dir,file)
            elif "_R2_" in file:
                pair2_file = os.path.join(sample_data_dir,file)
            else:
                sys.exit("Error: file {} does not respect the paired end format".format(file))

        outputFile1P = "{}_1.fastq.gz".format(sample_dir_name)
        outputFile1U = "{}_1_unp.fastq.gz".format(sample_dir_name)
        outputFile2P = "{}_2.fastq.gz".format(sample_dir_name)
        outputFile2U = "{}_2_unp.fastq.gz".format(sample_dir_name)

        currentCommand.append(pair1_file)
        currentCommand.append(pair2_file)
        currentCommand.append(outputFile1P)
        currentCommand.append(outputFile1U)
        currentCommand.append(outputFile2P)
        currentCommand.append(outputFile2U)

        if args.hard_trim == True:
            if args.keep == None:
                sys.exit("Error: if --hard-trim option is specified --keep must be spefified as well")
            currentCommand.append("CROP:{}".format(args.keep))
        currentCommand.append("LEADING:3")
        if args.no_quality_trim != True:
            currentCommand.append("SLIDINGWINDOW:4:15")
        currentCommand.append("MINLEN:{}".format(args.min_length))

        trimmomatic_stdOut = open("trimmomatic.stdOut", "a")
        trimmomatic_stdErr = open("trimmomatic.stdErr", "a")
        subprocess.call(currentCommand, stdout=trimmomatic_stdOut, stderr=trimmomatic_stdErr)


        command_singleReads = "zcat *unp.fastq.gz | gzip -c > {}_single.fastq.gz".format(sample_dir_name)
        subprocess.call(command_singleReads, stdout=trimmomatic_stdOut, stderr=trimmomatic_stdErr, shell=True)

        command_rm = ["rm", outputFile1U, outputFile2U]
        subprocess.call(command_rm, stdout=trimmomatic_stdOut, stderr=trimmomatic_stdErr)

        trimmomatic_stdOut.close()
        trimmomatic_stdErr.close()
        os.chdir(workingDir)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform hard trimming and quality trimming. Best suited for high numebr of samples. By default it removes leading \
        bases with low coverage and perform quality trimming. Output is returned in gz format')
    parser.add_argument('--sample-data-dir', help="full path to directory containing one folder per sample. Each sample contains only one library (i.e., one PE lib)", type=str)
    parser.add_argument('--hard-trim', help="perform only hard trim (--keep1 and keep2 must be specified)", action='store_true')
    parser.add_argument('--keep', help="keeps only the first keep bases in read 1 and read 2. (defult is no hard trim).", type=int)
    parser.add_argument('--no-quality-trim', help="no quality trim is performed", action='store_true')
    parser.add_argument('--min-length', help="minimum length to output", default=36,  type=int)
    parser.add_argument('--global-config', help="global configuration file with path to programs and options")
    parser.add_argument('--threads', help="Number of threads to use (overwrites global-config option) --> NOT WORKING")
    args = parser.parse_args()

    main(args)