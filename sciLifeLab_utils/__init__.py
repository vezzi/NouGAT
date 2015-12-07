from __future__ import absolute_import
from __future__ import print_function
import os
import yaml
import subprocess


def submit_job(sample_config, jobname, rundir, cliargs, extramodules=[]):
    """ Write a slurm file and sbatch it (if not dry-run)."""

    slurmfile_path = os.path.join(rundir, "{}.slurm".format(jobname))
    with open(slurmfile_path, "w") as slurmfile:
        slurmfile.write("#! /bin/bash -l\n")
        slurmfile.write("#SBATCH -A {}\n".format(cliargs.project))
        slurmfile.write("#SBATCH -o {}.out\n".format(jobname))
        slurmfile.write("#SBATCH -e {}.err\n".format(jobname))
        slurmfile.write("#SBATCH -J {}.job\n".format(jobname))
        if cliargs.threads<16 :
            slurmfile.write("#SBATCH -p core -n {}\n".format(cliargs.threads))
        else:
            slurmfile.write("#SBATCH -p node -n {}\n".format(cliargs.threads))
        slurmfile.write("#SBATCH -t {}\n".format(cliargs.time))
        if hasattr(cliargs, "email"):
            slurmfile.write("#SBATCH --mail-user {}\n".format(cliargs.email))
            slurmfile.write("#SBATCH --mail-type=ALL\n")
        if hasattr(cliargs, "qos"):
            slurmfile.write("#SBATCH --qos={}".format(cliargs.qos))
        slurmfile.write("\n\n")
        slurmfile.write("set -e\n")
        slurmfile.write("source activate {}\n".format(cliargs.env))
        slurmfile.write("module load bioinfo-tools\n")
        for module in extramodules:
            slurmfile.write(module)

        slurmfile.write("deNovo_pipeline.py --global-config {} "
            "--sample-config {}\n\n".format(cliargs.global_config, sample_config))

    command=("sbatch", slurmfile_path)
    print(command)
    try:
        if cliargs.dry_run:
            return 0
    except AttributeError as e:
        print("Warning! Could not determine if dry-run, running the command anyway: {}".format(e))
    return subprocess.call(command)
