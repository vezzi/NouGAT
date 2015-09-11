import sys, os, yaml, glob
import subprocess
import argparse
import re
from sciLifeLab_utils import submit_job

def main(args):
    projectFolder    = os.getcwd()
    samples_data_dir = args.sample_data_dir
    #UPPMAX assumption
    projectName      = os.path.basename(os.path.normpath(samples_data_dir))
    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) \
            if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        sample_folder = os.path.join(os.getcwd(), sample_dir_name)
        if not os.path.exists(sample_folder):
            os.makedirs(sample_folder)
        os.chdir(sample_folder)
        # if this is the case I need to retrive the project name from the yaml
        # file
        if args.afterqc:
            QC_YAML_file = os.path.join(samples_data_dir,sample_dir_name,
                    "{}_QCcontrol.yaml".format(sample_dir_name))
            if not os.path.exists(QC_YAML_file):
                sys.exit("Error file {} must exists!".format(QC_YAML_file))
            with open(QC_YAML_file) as QC_YAML_file_handle:
                QC_sample_config = yaml.load(QC_YAML_file_handle)
            # TODO: I need to use the sample sheet that must be present in
            # the QC folder to extract the project name
            projectName  = QC_sample_config["projectName"]
        #Now all the info is in place and I am in the correct folder
        pipeline = "assemble"
        tools = list(args.assemblers)
        tools = map(str, tools) # Beware whoever inputs unicode characters
        sample_YAML_name = os.path.join(sample_folder,
                "{}_{}.yaml".format(sample_dir_name, pipeline))
        sample_YAML = open(sample_YAML_name, 'w')
        sample_YAML.write("pipeline:\n")
        sample_YAML.write(" {}\n".format(pipeline))
        sample_YAML.write("tools:\n")
        sample_YAML.write(" {}\n".format(tools))
        sample_YAML.write("output: {}\n".format(sample_dir_name))
        sample_YAML.write("projectName: {}\n".format(projectName))
        sample_YAML.write("kmer: {}\n".format(args.kmer))
        sample_YAML.write("threads: {}\n".format(args.threads))
        sample_YAML.write("genomeSize: {}\n".format(args.genomesize))
        #I have to distinguish between afterQC and not
        sample_data_dir = ""
        sample_files = []
        if args.afterqc:
            sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
            fastq_files     = os.path.join(sample_data_dir, "results",
                    "fastq_trimmed")
            sample_files    = [os.path.join(fastq_files, f) for f in \
                    os.listdir(fastq_files) \
                    if (os.path.isfile(os.path.join(fastq_files,f)) \
                    and re.search('[1|2].fastq.gz$',f))]
        else:
            sample_data_dir = os.path.join(samples_data_dir,sample_dir_name)
            # full path to flowcell
            flowcells_dirs  = [os.path.join(sample_data_dir,flowcell) \
                    for flowcell in os.listdir(sample_data_dir) \
                    if os.path.isdir(os.path.join(sample_data_dir,flowcell))]
            for flowcell in flowcells_dirs:
                sample_files.extend([os.path.join(flowcell, f) for f in \
                        os.listdir(flowcell) \
                        if (os.path.isfile(os.path.join(flowcell,f)) \
                        and re.search('.gz$',f))])
        # now sample_files contains all the file sequenced for this sample
        pair1_file = ""
        pair2_file = ""
        single = ""
        library = 1
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

        #Run the job
        all_modules = {
                "abyss": "module load abyss/1.3.5\n",
                "soapdenovo": "module load soapdenovo/2.04-r240\n",
                "spades": "module load spades/3.6.0\n",
                "cabog": "module load cabog/8.1\n",
                "allpaths": "module unload gcc\nmodule load allpathslg/52485\n"
        }

        extramodules = [all_modules[tool] for tool in tools]
        jobname = "{}_{}".format(sample_dir_name, pipeline)
        submit_job(sample_YAML_name, jobname, os.getcwd(), args, extramodules)
        os.chdir(projectFolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This is a utility script to '
            'generate de novo assemblies with sample sequenced with a single '
            'library. For samples sequenced with multiple libraries the user '
            'needs to prepare the sample_config file and run the '
            'deNovoPipeline manually. If multiple samples are present in the '
            'specified folder then one assembly for each sample will be '
            'performed. If a sample is splitted across multiple runs all the '
            'data willl be used')
    parser.add_argument('--global-config', type=str, required=True,
            help="global configuration file")
    parser.add_argument('--sample-data-dir', type=str, required=True,
            help=("Path to directory (usually INBOX or QC output) containing "
            "the project (in INBOX case scilife structure "
            "project/sample/flowcell/ is assumed)"))
    parser.add_argument('--orientation', type=str, required=True,
            help="orientation of the libraries")
    parser.add_argument('--insert', type=str, required=True,
            help="expected insert size of the libraries")
    parser.add_argument('--std', type=str, required=True,
            help=("expected stdandard variation of the insert size of the "
            "libraries"))
    parser.add_argument('--env', type=str,
            default="DeNovoPipeline", help=("name of the virtual enviorment "
            "(default is DeNovoPipeline)"))
    parser.add_argument('--assemblers', type=str, required=True,
            nargs='+', help=("List of assemblers to be "
            "employed on the datasets specified"))
    parser.add_argument('--kmer', type=int, required=True,
            help=("kmer size to employ when requested (i.e., some tools "
            "make a guess"))
    parser.add_argument('--genomesize', type=int, required=True,
            help="Estimated genome size (make an educated guess)")
    parser.add_argument('--afterqc', action='store_true',
            default = False,  help=("To be specified if sample-data-dir is a "
            "QC output"))
    parser.add_argument('--email', type=str,
            help="Send notifications/job status updates to this email address")
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




