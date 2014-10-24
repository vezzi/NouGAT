import sys, os, yaml, glob
import subprocess
import argparse
import re


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
        if args.reference == "":
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

        if args.reference is not "":
            sample_YAML.write("reference: {}\n".format(args.reference))
        sample_YAML.write("libraries:\n")

        sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
        # full path to flowcell
        flowcells_dirs  = [os.path.join(sample_data_dir,flowcell) \
                for flowcell in os.listdir(sample_data_dir) \
                if os.path.isdir(os.path.join(sample_data_dir,flowcell))]

        sample_files = []
        for flowcell in flowcells_dirs:
            sample_files.extend([os.path.join(flowcell, f) for f in \
                    os.listdir(flowcell) \
                    if (os.path.isfile(os.path.join(flowcell,f)) \
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
        if not hasattr(args, "email"):
            args.email = None
        submit_job(sample_YAML_name, args.global_config, sample_dir_name ,
                pipeline, args.env, args.email, args.time, args.project,
                args.threads, args.qos) # now I can submit the job to slurm
        os.chdir(projectFolder)


def submit_job(sample_config, global_config, output,  pipeline, env,
        email=None, required_time='1-00:00:00', project='a2010002',
        threads=16, qos=None):

    workingDir = os.getcwd()
    slurm_file = os.path.join(workingDir, "{}_{}.slurm".format(
        output,pipeline))
    slurm_handle = open(slurm_file, "w")

    slurm_handle.write("#! /bin/bash -l\n")
    slurm_handle.write("set -e\n")
    slurm_handle.write("#SBATCH -A {}\n".format(project))
    slurm_handle.write("#SBATCH -o {}_{}.out\n".format(output,pipeline))
    slurm_handle.write("#SBATCH -e {}_{}.err\n".format(output,pipeline))
    slurm_handle.write("#SBATCH -J {}_{}.job\n".format(output,pipeline))
    if threads<16 :
        slurm_handle.write("#SBATCH -p core -n {}\n".format(threads))
    else:
        slurm_handle.write("#SBATCH -p node -n {}\n".format(threads))
    slurm_handle.write("#SBATCH -t {}\n".format(required_time))
    if email:
        slurm_handle.write("#SBATCH --mail-user {}\n".format(email))
        slurm_handle.write("#SBATCH --mail-type=ALL\n")
    if qos:
        slurm_handle.write("#SBATCH --qos={}".format(qos))
    slurm_handle.write("\n\n")
    slurm_handle.write("source activate {}\n".format(env))
    slurm_handle.write("module load bioinfo-tools\n")
    slurm_handle.write("module load bwa\n")
    slurm_handle.write("module load abyss/1.3.5\n")
    slurm_handle.write("deNovo_pipeline.py --global-config {} "
            "--sample-config {}\n\n".format(global_config,sample_config))
    slurm_handle.close()

    command=("sbatch", slurm_file)
    print command
    subprocess.call(command) 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=str, default="",
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
    parser.add_argument('--qos', type=str, default=None,
            help=("Specify a quality of service preset for the job (eg. "
            "--qos short)"))
    args = parser.parse_args()

    main(args)
