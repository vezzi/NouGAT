import sys, os, yaml, glob
import subprocess
import argparse
import re


def main(args):
    projectFolder    = os.getcwd()
    samples_data_dir = args.sample_data_dir
    projectName      = os.path.basename(os.path.normpath(samples_data_dir)) #UPPMAX assumption
    
    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        sample_folder = os.path.join(os.getcwd(), sample_dir_name)
        if not os.path.exists(sample_folder):
            os.makedirs(sample_folder)
        os.chdir(sample_folder)
        if args.afterQC: #if this is the case I need to retrive the project name from the yaml file
            QC_YAML_file = os.path.join(samples_data_dir,sample_dir_name, "{}_QCcontrol.yaml".format(sample_dir_name))
            if not os.path.exists(QC_YAML_file):
                sys.exit("Error file {} must exists!".format(QC_YAML_file))
            with open(QC_YAML_file) as QC_YAML_file_handle:
                QC_sample_config = yaml.load(QC_YAML_file_handle)
            projectName  = QC_sample_config["projectName"] # TODO: I need to use the sample sheet that must be present in the QC folder to extract the project name
        #Now all the info is in place and I am in the correct folder
        pipeline = "assemble"
        tools    = args.assemblers
        sample_YAML_name = os.path.join(sample_folder,  "{}_{}.yaml".format(sample_dir_name, pipeline))
        sample_YAML      = open(sample_YAML_name, 'w')
        sample_YAML.write("pipeline:\n")
        sample_YAML.write(" {}\n".format(pipeline))
        sample_YAML.write("tools:\n")
        sample_YAML.write(" {}\n".format(tools[0]))
        sample_YAML.write("output: {}\n".format(sample_dir_name))
        sample_YAML.write("projectName: {}\n".format(projectName))
        sample_YAML.write("kmer: {}\n".format(args.kmer))
        sample_YAML.write("threads: {}\n".format(args.threads))
        sample_YAML.write("genomeSize: {}\n".format(args.genomeSize))
        #I have to distinguish between afterQC and not
        sample_data_dir = ""
        sample_files = []
        if args.afterQC:
            sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
            fastq_files     = os.path.join(sample_data_dir, "results", "fastq_trimmed")
            sample_files    = [os.path.join(fastq_files, f) for f in os.listdir(fastq_files) if (os.path.isfile(os.path.join(fastq_files,f)) and re.search('[1|2].fastq.gz$',f))]
        else:
            sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
            flowcells_dirs  = [os.path.join(sample_data_dir,flowcell)  for flowcell in os.listdir(sample_data_dir) if os.path.isdir(os.path.join(sample_data_dir,flowcell))] # full path to flowcell
            for flowcell in flowcells_dirs:
                sample_files.extend([os.path.join(flowcell, f) for f in os.listdir(flowcell) if (os.path.isfile(os.path.join(flowcell,f)) and re.search('.gz$',f))])
        # now sample_files contains all the file sequenced for this sample
        pair1_file = ""
        pair2_file = ""
        single     = ""
        library    = 1
        sample_YAML.write("libraries:\n")
        for file in sample_files:
            sample_YAML.write(" lib{}:\n".format(library))
            if "_1.fastq.gz" in file:
                pair1_file = file
                pair2_file = re.sub("_1.fastq.gz", "_2.fastq.gz", file)
            elif "_2.fastq.gz" in file:
                pair2_file = file
                pair1_file = re.sub("_2.fastq.gz", "_1.fastq.gz", file)

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
        submit_job(sample_YAML_name, args.global_config, sample_dir_name , pipeline, args.env, args.email, args.time, args.project, args.threads) # now I can submit the job to slurm
#        submit_job(sample_YAML_name, args.global_config, sample_dir_name , pipeline, args.env) # now I can submit the job to slurm
        os.chdir(projectFolder)

def submit_job(sample_config, global_config, output,  pipeline, env, email=None, required_time='1-00:00:00', project='a2010002', threads=16 ):
    workingDir = os.getcwd()
    slurm_file = os.path.join(workingDir, "{}_{}.slurm".format(output,pipeline))
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
    slurm_handle.write("\n\n")
    slurm_handle.write("source activate {}\n".format(env))
    slurm_handle.write("module load bioinfo-tools\n")
    slurm_handle.write("module load soapdenovo/2.04-r240\n")
    slurm_handle.write("module load spades/3.0.0\n")
    slurm_handle.write("module load cabog/8.1\n")
    slurm_handle.write("module load abyss/1.3.5\n")
    slurm_handle.write("deNovo_pipeline.py --global-config {} --sample-config {}\n\n".format(global_config,sample_config))
    slurm_handle.close()
    command=("sbatch", slurm_file)
    print command
    subprocess.call(command)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a utility script to generate de novo assemblies with sample sequenced with a single library. For samples sequenced\
with multiple libraries the user needs to prepare the sample_config file and run the deNovoPipeline manually. If multiple samples are present in the specified folder then one\
assembly for each sample will be performed. If a sample is splitted across multiple runs all the data willl be used')
    parser.add_argument('--global-config'  , type=str, required=True, help="global configuration file")
    parser.add_argument('--sample-data-dir', type=str, required=True, help="Path to directory (usually INBOX or QC output) containing the project (in INBOX case scilife structure project/sample/flowcell/ is assumed)")
    parser.add_argument('--orientation'    , type=str, required=True, help="orientation of the libraries")
    parser.add_argument('--insert'         , type=str, required=True, help="expected insert size of the libraries")
    parser.add_argument('--std'            , type=str, required=True, help="expected stdandard variation of the insert size of the libraries")
    parser.add_argument('--env'            , type=str, default="DeNovoPipeline", help="name of the virtual enviorment (default is DeNovoPipeline)")
    parser.add_argument('--assemblers'     , type=str, required=True, action='append', nargs='+', help="List of assemblers to be employed on the datasets specified")
    parser.add_argument('--kmer'           , type=int, required=True, help="kmer size to employ when requested (i.e., some tools make a guess")
    parser.add_argument('--genomeSize'     , type=int, required=True,  help="Estimated genome size (make an educated guess)")
    parser.add_argument('--afterQC'        , action='store_true', default = False,  help="To be specified if sample-data-dir is a QC output")
    parser.add_argument('--email'          , type=str, help="Send notifications/job status updates to this email address.")
    parser.add_argument('--time'           , type=str, default="1-00:00:00", help="required time for the job (default is 1 day : 1-00:00:00)")
    parser.add_argument('--project'        , type=str, default="a2010002", help="project name for slurm submission (default is a2010002)")
    parser.add_argument('--threads'        , type=int, default=16, help="Number of thread the job will require")
    args = parser.parse_args()
    main(args)




