import sys, os, yaml, glob
import subprocess
import argparse



def main(args):
    workingDir = os.getcwd()
    samples_data_dir  = args.sample_data_dir
    assemblies_data_dir = args.assemblies_data_dir
    assemblers = sum(args.assemblers, [])
    
    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        validation_folder = os.path.join(os.getcwd(), sample_dir_name)
        if not os.path.exists(validation_folder):
            os.makedirs(validation_folder)
        else:
            print "done ({} folder already present, assumed already run)".format(validation_folder)
        os.chdir(validation_folder)
    
        sample_YAML_name = os.path.join(validation_folder, "{}_validation.yaml".format(sample_dir_name))
        sample_YAML = open(sample_YAML_name, 'w')
       
        sample_YAML.write("pipeline:\n")
        sample_YAML.write(" evaluete\n")
        sample_YAML.write("tools:\n")
        sample_YAML.write(" [qaTools, FRC]\n")
        sample_YAML.write("genomeSize: {}\n".format(args.genomeSize))
        sample_YAML.write("output: {}\n".format(sample_dir_name))
        sample_YAML.write("minCtgLength: 2000\n")
        sample_YAML.write("reference:\n")
        sample_YAML.write("libraries:\n")
        
        sample_data_dir   = os.path.join(samples_data_dir,sample_dir_name)
        assembly_data_dir = os.path.join(assemblies_data_dir,sample_dir_name)
        sample_files = [ f for f in os.listdir(sample_data_dir) if os.path.isfile(os.path.join(sample_data_dir,f))]

        pair1_file = ""
        pair2_file = ""
        single     = ""
        sample_YAML.write(" lib1:\n")
        for file in sample_files:
            if "_R1_" in file or "_1.fastq.gz" in file:
                if pair1_file:
                    sys.exit("Error: processing sample {} found more that one library/run for read 1".format(sample_dir_name))
                pair1_file = os.path.join(sample_data_dir,file)
                sample_YAML.write("  pair1: {}\n".format(pair1_file))
            elif "_R2_" in file or "_2.fastq.gz" in file:
                if pair2_file:
                    sys.exit("Error: processing sample {} found more that one library/run for read 2".format(sample_dir_name))
                pair2_file = os.path.join(sample_data_dir,file)
                sample_YAML.write("  pair2: {}\n".format(pair2_file))
            elif "merged" in file or "single" in file:
                single = os.path.join(sample_data_dir,file)

        sample_YAML.write("  orientation: {}\n".format(args.orientation))
        sample_YAML.write("  insert: {}\n".format(args.insert))
        sample_YAML.write("  std: {}\n".format(args.std))

        sample_YAML.close()
        
        command_to_run = "python ~/assembly_pipeline/de_novo_scilife/utils/run_assembly_evaluation.py --global-config {} --sample-config {}_validation.yaml \
 --assemblies-dir {} --assembler {}".format(args.global_config, sample_dir_name, assembly_data_dir , " ".join(assemblers))
        print command_to_run
        subprocess.call(command_to_run, shell=True)
        
        
        os.chdir(workingDir)


def submit_job(sample_config, global_config, sample_name):
    workingDir = os.getcwd()
    slurm_file = os.path.join(workingDir, "{}.slurm".format(sample_name))
    slurm_handle = open(slurm_file, "w")

    slurm_handle.write("#! /bin/bash -l\n")
    slurm_handle.write("set -e\n")
    slurm_handle.write("#SBATCH -A a2010002\n")
    slurm_handle.write("#SBATCH -o {}_QC.out\n".format(sample_name))
    slurm_handle.write("#SBATCH -e {}_QC.err\n".format(sample_name))
    slurm_handle.write("#SBATCH -J {}_QC.job\n".format(sample_name))
    slurm_handle.write("#SBATCH -p node -n 8\n")
    slurm_handle.write("#SBATCH -t 05:00:00\n")
    slurm_handle.write("#SBATCH --mail-user francesco.vezzi@scilifelab.se\n")
    slurm_handle.write("#SBATCH --mail-type=ALL\n")
    
    slurm_handle.write("\n\n");
    slurm_handle.write("module load abyss/1.3.5\n");
    slurm_handle.write("python ~/assembly_pipeline/de_novo_scilife/script/deNovo_pipeline.py --global-config {} --sample-config {}\n\n".format(global_config,sample_config))
    slurm_handle.close()
    
    command=("sbatch", slurm_file)
    print command
    subprocess.call(command)
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-data-dir', help="full path to directory containing one folder per sample. Each sample contaains only one library (i.e., one PE lib)", type=str)
    parser.add_argument('--assemblies-data-dir', help="full path to directory containing one for each sample the de novo assemblies", type=str)
    parser.add_argument('--assemblers', action='append', nargs='+', help="List of assemblers to be evalueted")
    parser.add_argument('--genomeSize', help="genome size", type=str)
    parser.add_argument('--orientation', help="I assume I am working only with PE (if not manual editing is needed)", type=str)
    parser.add_argument('--insert', help="I assume that all samples have the same insert (if not manual editing is needed)", type=str)
    parser.add_argument('--std', help="I assume tha all sample have same std (if not manual editing is needed)", type=str)
    parser.add_argument('--global-config', help="global configuration file")
    
    args = parser.parse_args()

    main(args)