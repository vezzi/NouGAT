import argparse
import os
import glob
import re
import subprocess

def main(arg):
    project = os.path.split(os.path.realpath(arg.source))[1]
    move_from_path = "{}/*/results/".format(arg.source)
    pathway = glob.glob(move_from_path)
    samples = []
    for path in pathway:
        try:
            samples.append(os.path.basename( os.path.dirname(path.rstrip('/'))))

        except IndexError:
            pass
    
    
    for sample in samples:
        dest = "/proj/{}/INBOX/{}/QC_analysis/{}".format(arg.uppnexid,project,sample)
        try:
             os.mkdir("/proj/{}/INBOX/{}/QC_analysis".format(arg.uppnexid,project))
          
        except OSError:
            pass
        try:
           os.mkdir(dest)
       
        except OSError:
            pass
       
        cmd = ["rsync", "-auhvr", "{}/{}/results".format(arg.source, sample), dest]
        subprocess.call(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type = str,
            help= "Path to QC analysis folder")
    parser.add_argument("--uppnexid", required=True, type = str,
            help =("Destination Uppnex id"))
    projectID = parser.parse_args()
    main(projectID)


