import sys, os, yaml, glob
import subprocess
import argparse
import re
from sciLifeLab_utils import submit_job

def main(args):

    projectFolder = os.getcwd()
    assemblies_dir  = args.assembly_dir

    for sample_dir_name in [dir for dir in os.listdir(assemblies_dir) \
            if os.path.isdir(os.path.join(assemblies_dir, dir))]:
        # in this folder I stored all the assemblies
        assemblies_folder   = os.path.join(assemblies_dir, sample_dir_name)
        # in this folder I will compute the validation
        validation_folder = os.path.join(os.getcwd(), sample_dir_name)
        if not os.path.exists(validation_folder):
            os.makedirs(validation_folder)
        os.chdir(validation_folder)
        ## Restore all information present in the sample yaml file present
        ## in the assembly folder
        sample_config_assembly_file = os.path.join(assemblies_folder,
                "{}_assemble.yaml".format(sample_dir_name))
        with open(sample_config_assembly_file) as sample_config_handle:
            sample_config_assembly = yaml.load(sample_config_handle)
        # prepare for each assembly employed a validation job --- do not check
        # in sample sheet the asse,blies, run one validation per forder
        # only (assumptions are done on assmebly name)
        for assembler in [dir for dir in os.listdir(assemblies_folder) \
                if os.path.isdir(os.path.join(assemblies_folder, dir))]:
            if not os.path.exists(assembler):
                os.makedirs(assembler)
            os.chdir(assembler)
            assembly_dir = os.path.join(assemblies_folder, assembler)
            assembly_name = os.path.join(assembly_dir,
                    "{}.scf.fasta".format(sample_config_assembly["output"]))
            pipeline = "evaluete"
            sample_YAML_name = "{}_{}.yaml".format(sample_dir_name, pipeline)
            sample_YAML = open(sample_YAML_name, 'w')
            sample_YAML.write("pipeline:\n")
            sample_YAML.write(" {}\n".format(pipeline))
            sample_YAML.write("tools:\n")
            sample_YAML.write(" [align, qaTools, FRC]\n")
            sample_YAML.write(
                    "output: {}\n".format(sample_config_assembly["output"]))
            sample_YAML.write(
                    "projectName: {}_validate\n".format(
                    sample_config_assembly["output"]))
            sample_YAML.write("kmer: \n")
            sample_YAML.write("threads: {}\n".format(args.threads))
            sample_YAML.write(
                    "genomeSize: {}\n".format(
                    sample_config_assembly["genomeSize"]))
            sample_YAML.write("minCtgLength: 1000\n")
            sample_YAML.write("reference: {}\n".format(assembly_name))
            sample_YAML.write("libraries:\n")
            for library, libraryData in \
                    sample_config_assembly["libraries"].items():
                sample_YAML.write(" {}:\n".format(library))
                sample_YAML.write("  pair1: {}\n".format(libraryData["pair1"]))
                sample_YAML.write("  pair2: {}\n".format(libraryData["pair2"]))
                sample_YAML.write("  orientation: {}\n".format(
                    libraryData["orientation"]))
                sample_YAML.write("  insert: {}\n".format(
                    libraryData["insert"]))
                sample_YAML.write("  std: {}\n".format(libraryData["std"]))

            sample_YAML.close

            # Run the job
            extramodules = ["module load samtools/1.1\nmodule load bwa\n"]
            jobname = "{}_{}_{}".format(sample_dir_name, pipeline, assembler)
            submit_job(sample_YAML_name, jobname, os.getcwd(), args, extramodules)

            os.chdir(validation_folder)
        os.chdir(projectFolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is a utility script to "
            "run assembly validation with sample sequenced with a single "
            "library. For samples sequenced with multiple libraries the user "
            "needs to prepare the sample_config file and run the "
            "deNovoPipeline manually. If multiple samples are present in the "
            "specified folder then one assembly for each sample will be "
            "performed. If a sample is splitted across multiple runs all the "
            "data willl be used")
    parser.add_argument('--global-config', type=str, required=True,
            help="global configuration file")
    parser.add_argument('--assembly-dir', type=str, required=True,
            help="Path to directory containg assemblies. All meta-deta is "
            "extracted from sample config file present in each sample folder "
            "(this behaviour is over-written by --local-config)")
    parser.add_argument('--env', type=str, default="DeNovoPipeline",
            help="name of the virtual enviorment (default is DeNovoPipeline)")
    parser.add_argument('--email', type=str, default=None,
            help="Send notifications/job status updates to this email "
            "address.")
    parser.add_argument('--time', type=str, default="1-00:00:00",
            help="required time for the job (default is 1 day : 1-00:00:00)")
    parser.add_argument('--project', type=str, default="a2010002",
            help="project name for slurm submission (default is a2010002)")
    parser.add_argument('--threads', type=int, default=16,
            help="Number of thread the job will require")
    parser.add_argument('--multiple-lib-proj', action='store_true',
            default = False,  help="To be specified if we are running a "
            "mulitple library assembly")
    parser.add_argument('--qos', type=str, default=None,
            help="Specify a quality of service preset for the job "
            "(eg. --qos short)")
    args = parser.parse_args()
    main(args)




