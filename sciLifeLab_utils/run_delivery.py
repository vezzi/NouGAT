import argparse
import os
import glob
import re
import subprocess

def main(arg):
    project = os.path.split(os.path.realpath(arg.source))[1]
    move_from_path = "{}/*/".format(arg.source)
    samples_path = glob.glob(move_from_path)
    #create in inbox folder for project
    #dest = "/proj/{}/INBOX/{}/".format(arg.uppnexid,project)
    dest = "/Users/vezzi/Desktop/"
    try:
        os.makedirs(dest)
    except OSError:
        pass
    for sample in samples_path:
        sampel_name = os.path.basename(os.path.normpath(sample))
        sample_dest = os.path.join(dest, sampel_name)
        cmd = ["rsync", "-auhvr",  sample, sample_dest]
        subprocess.call(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type = str,
            help= "Path to QC analysis folder")
    parser.add_argument("--uppnexid", required=True, type = str,
            help =("Destination Uppnex id"))
    projectID = parser.parse_args()
    main(projectID)


