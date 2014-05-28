import sys, os, yaml, glob
import subprocess
import argparse
import re


def main(args):
    projectFolder = os.getcwd()
    assemblies_dir  = args.assembly_dir
    for sample_dir_name in [dir for dir in os.listdir(assemblies_dir) if os.path.isdir(os.path.join(assemblies_dir, dir))]:
        assemblies_folder   = os.path.join(assemblies_dir, sample_dir_name) # in this folder I stored all the assemblies
        validation_folder = os.path.join(os.getcwd(), sample_dir_name)  # in this folder I will compute the validation
        if not os.path.exists(validation_folder):
            os.makedirs(validation_folder)
        os.chdir(validation_folder)
        ## Restore all information present in the sample yaml file present in the assembly folder
        sample_config_assembly_file = os.path.join(assemblies_folder, "{}_assemble.yaml".format(sample_dir_name))
        with open(sample_config_assembly_file) as sample_config_handle:
            sample_config_assembly = yaml.load(sample_config_handle)
        #prepare for each assembly employed a validation job --- do not check in sample sheet the asse,blies, run one validation per forder only (assumptions are done on assmebly name)
        for assembler in [dir for dir in os.listdir(assemblies_folder) if os.path.isdir(os.path.join(assemblies_folder, dir))]:
            if not os.path.exists(assembler):
                os.makedirs(assembler)
            os.chdir(assembler)
            assembly_dir  = os.path.join(assemblies_folder, assembler)
            assembly_name = os.path.join(assembly_dir, "{}.scf.fasta".format(sample_config_assembly["output"]))
            pipeline = "evaluete"
            sample_YAML_name = "{}_{}.yaml".format(sample_dir_name, pipeline)
            sample_YAML      = open(sample_YAML_name, 'w')
            sample_YAML.write("pipeline:\n")
            sample_YAML.write(" {}\n".format(pipeline))
            sample_YAML.write("tools:\n")
            sample_YAML.write(" [align, qaTools, FRC]\n")
            sample_YAML.write("output: {}\n".format(sample_config_assembly["output"]))
            sample_YAML.write("projectName: {}_validate\n".format(sample_config_assembly["output"]))
            sample_YAML.write("kmer: \n")
            sample_YAML.write("threads: 16\n")
            sample_YAML.write("genomeSize: {}\n".format(sample_config_assembly["genomeSize"]))
            sample_YAML.write("minCrgLength: 1000\n")
            sample_YAML.write("reference: {}\n".format(assembly_name))
            sample_YAML.write("libraries:\n")
            for library, libraryData in sample_config_assembly["libraries"].items():
                sample_YAML.write(" {}:\n".format(library))
                sample_YAML.write("  pair1: {}\n".format(libraryData["pair1"]))
                sample_YAML.write("  pair2: {}\n".format(libraryData["pair2"]))
                sample_YAML.write("  orientation: {}\n".format(libraryData["orientation"]))
                sample_YAML.write("  insert: {}\n".format(libraryData["insert"]))
                sample_YAML.write("  std: {}\n".format(libraryData["std"]))

            sample_YAML.close
            if not hasattr(args, "email"):
                args.email=None

            submit_job(sample_YAML_name, args.global_config, sample_dir_name , pipeline, assembler, args.env, args.email) # now I can submit the job to slurm
            os.chdir(validation_folder)
        os.chdir(projectFolder)

def submit_job(sample_config, global_config, output,  pipeline, assembler, env, email=None):
    workingDir = os.getcwd()
    slurm_file = os.path.join(workingDir, "{}_{}_{}.slurm".format(output,pipeline, assembler))
    slurm_handle = open(slurm_file, "w")
    slurm_handle.write("#! /bin/bash -l\n")
    slurm_handle.write("set -e\n")
    slurm_handle.write("#SBATCH -A b2013064\n")
    slurm_handle.write("#SBATCH -o {}_{}_{}.out\n".format(output,pipeline,assembler))
    slurm_handle.write("#SBATCH -e {}_{}_{}.err\n".format(output,pipeline,assembler))
    slurm_handle.write("#SBATCH -J {}_{}_{}.job\n".format(output,pipeline,assembler))
    slurm_handle.write("#SBATCH -p node -n 16\n")
    slurm_handle.write("#SBATCH -t 1-00:00:00\n")
    if email :
        slurm_handle.write("#SBATCH --mail-user {}\n".format(email))
        slurm_handle.write("#SBATCH --mail-type=ALL\n")

    slurm_handle.write("\n\n")
    slurm_handle.write("source activate {}\n".format(env))
    slurm_handle.write("load_modules\n")
    slurm_handle.write("module load abyss/1.3.5\n")
    slurm_handle.write("deNovo_pipeline.py --global-config {} --sample-config {}\n\n".format(global_config,sample_config))
    slurm_handle.close()
    command=("sbatch", slurm_file)
    print command
    subprocess.call(command)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a utility script to run assembly validation with sample sequenced with a single library. For samples sequenced\
with multiple libraries the user needs to prepare the sample_config file and run the deNovoPipeline manually. If multiple samples are present in the specified folder then one\
assembly for each sample will be performed. If a sample is splitted across multiple runs all the data willl be used')
    parser.add_argument('--global-config'  , type=str, required=True, help="global configuration file")
    parser.add_argument('--assembly-dir'   , type=str, required=True, help="Path to directory containg assemblies. All meta-deta is extracted from sample config file present in each sample folder")
    parser.add_argument('--env'            , type=str, default="DeNovoPipeline", help="name of the virtual enviorment (default is DeNovoPipeline)")
    parser.add_argument('--email'	   , type=str, default=None, help="Send notifications/job status updates to this email address.")
    args = parser.parse_args()
    main(args)




