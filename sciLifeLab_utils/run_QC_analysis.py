from __future__ import absolute_import
import sys, os, yaml, glob
import subprocess
import argparse
import re
from sciLifeLab_utils import submit_job

def main(args):
    projectFolder    = os.getcwd()
    samples_data_dir = args.sample_data_dir
    projectName      = os.path.basename(os.path.normpath(samples_data_dir))
    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) \
            if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        sample_folder = os.path.join(os.getcwd(), sample_dir_name)
        if not os.path.exists(sample_folder):
            os.makedirs(sample_folder)
        os.chdir(sample_folder)
        # now I am in the folder, i can run at the same time QC and MP anlaysis

        pipeline = "QCcontrol"
        tools    = ["trimmomatic", "fastqc", "abyss", "align"]
        if args.reference is None:
            tools    = ["trimmomatic", "fastqc", "abyss"]

        sample_YAML_name = os.path.join(sample_folder,  "{}_{}.yaml".format(
            sample_dir_name, pipeline))
        sample_YAML = open(sample_YAML_name, 'w')

        sample_YAML.write("pipeline:\n")
        sample_YAML.write(" {}\n".format(pipeline))
        sample_YAML.write("tools:\n")
        sample_YAML.write(" {}\n".format(tools))
        ##TODO: output must became sampleName
        sample_YAML.write("output: {}\n".format(sample_dir_name))
        sample_YAML.write("projectName: {}\n".format(projectName))
        sample_YAML.write("kmer: 35\n")
        sample_YAML.write("threads: {}\n".format(args.threads))
        sample_YAML.write("genomeSize: \n")
        sample_YAML.write("adapters: {}\n".format(args.adapter))

        if args.reference is not None:
            sample_YAML.write("reference: {}\n".format(args.reference))
        sample_YAML.write("libraries:\n")

        sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
        # full path to flowcell
        flowcells_dirs  = [os.path.join(sample_data_dir,flowcell) \
                for flowcell in os.listdir(sample_data_dir) \
                if os.path.isdir(os.path.join(sample_data_dir,flowcell))]
        
        # to adapt the diretory structure in IRMA where it have one folder for lib prep.
        # So shcek if the collected FC directories are really FC directory or a LIB prep dir
        lib_prep_dirs = [fc for fc in flowcells_dirs if re.match(r'^[A-Z]$', os.path.basename(fc))]
        for prep_dir in lib_prep_dirs:
            flowcells_dirs.extend([os.path.join(prep_dir,flowcell) \
                for flowcell in os.listdir(prep_dir) \
                if os.path.isdir(os.path.join(prep_dir,flowcell))])
            if prep_dir in flowcells_dirs:
                flowcells_dirs.remove(prep_dir)
            
        sample_files = []
        for flowcell in flowcells_dirs:

            sample_files.extend([os.path.join(flowcell, f) for f in \
                    os.listdir(flowcell) \
                    if (os.path.isfile(os.path.realpath(os.path.join(flowcell,f))) \
                    and re.search('.gz$',f))])
        # now sample_files contains all the file sequenced for this sample
        pair1_file = ""
        pair2_file = ""
        single     = ""
        library    = 1
        while len(sample_files) > 0:
            file = sample_files[0]
            sample_YAML.write(" lib{}:\n".format(library))
            if "_1.fastq.gz" in file:
                pair1_file = file
                pair2_file = re.sub("_1.fastq.gz", "_2.fastq.gz", file)
            elif "_2.fastq.gz" in file:
                pair2_file = file
                pair1_file = re.sub("_2.fastq.gz", "_1.fastq.gz", file)
            elif "R1_001.fastq.gz" in file:
                pair1_file = file
                pair2_file = re.sub("R1_001.fastq.gz", "R2_001.fastq.gz", file)
            elif "R2_001.fastq.gz" in file:
                pair2_file = file
                pair1_file = re.sub("R2_001.fastq.gz", "R1_001.fastq.gz", file)
            else:
                sys.exit("file {} does not respect naming convection. \
                        Exit!".format(file))

            sample_YAML.write("  pair1: {}\n".format(pair1_file))
            sample_YAML.write("  pair2: {}\n".format(pair2_file))
            sample_YAML.write("  orientation: {}\n".format(args.orientation))
            sample_YAML.write("  insert: {}\n".format(args.insert))
            sample_YAML.write("  std: {}\n".format(args.std))
            sample_files.remove(pair1_file)
            sample_files.remove(pair2_file)
            library += 1
        sample_YAML.close

        # Run the job
        extramodules = []
        if "abyss" in tools:
            extramodules.append("module load abyss/1.3.5\n")
        if "align" in tools:
            extramodules.append("module load samtools\nmodule load bwa\n")
        jobname = "{}_{}".format(sample_dir_name, pipeline)
        submit_job(sample_YAML_name, jobname, os.getcwd(), args, extramodules)
        os.chdir(projectFolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=str, default=None,
            help="path to the reference file")
    parser.add_argument('--adapter', type=str, required=True,
            help="path to the file containing the adaptor sequence to be removed")
    parser.add_argument('--global-config', type=str, required=True,
            help="global configuration file")
    parser.add_argument('--sample-data-dir', type=str, required=True,
            help=("Path to directory (usually INBOX) containing the project "
            "(one dir per sample, scilife structure project/sample/flowcell/)"))
    parser.add_argument('--orientation', type=str, required=True,
            help="orientation of the libraries")
    parser.add_argument('--insert', type=str, required=True,
            help="expected insert size of the libraries")
    parser.add_argument('--std', type=str, required=True,
            help=("expected stdandard variation of the insert size of "
            "the libraries"))
    parser.add_argument('--env', type=str,
            default="DeNovoPipeline", help=("name of the virtual enviorment "
            "(default is DeNovoPipeline)"))
    parser.add_argument('--email', type=str, 
            help=("Send notifications/job status updates to this email "
            "address."))
    parser.add_argument('--time', type=str, default="1-00:00:00",
            help="required time for the job (default is 1 day : 1-00:00:00)")
    parser.add_argument('--project', type=str, default="a2010002",
            help="project name for slurm submission (default is a2010002)")
    parser.add_argument('--threads', type=int, default=16,
            help="Number of thread the job will require")
    parser.add_argument('--qos', type=str,
            help=("Specify a quality of service preset for the job (eg. "
            "--qos short)"))
    args = parser.parse_args()

    main(args)
