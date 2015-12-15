from __future__ import absolute_import
from __future__ import print_function
import sys, os, yaml, glob
import argparse
from nougat import evaluete, assemble, QCcontrol, align, common


def main(args):
    with open(args.global_config) as in_handle:
        global_config = yaml.load(in_handle)
    with open(args.sample_config) as sample_config_handle:
        sample_config = yaml.load(sample_config_handle)

    check_consistency(global_config,sample_config)
    if common.check_dryrun(sample_config):
        print("Option dryrun idenitfied: commands will only be printed,)",
               "not executed")

    if sample_config["pipeline"] in global_config["Pipelines"]:
        run_analys(global_config, sample_config)
    else:
        sys.exit("Error: pipeline {} is not one of the supported \
                ones:{}".format(sample_config["pipeline"], 
                global_config["Pipelines"]))

    return 0


def run_analys(global_config, sample_config):
    """ check that user specified commands are supported by this pipeline and 
    that all commands I am going to run are available either via PATH or 
    via global config"""
    # check that what I am going to run is available on the path and on the 
    # global config
    common._check_pipeline(sample_config, global_config)
    pipeline = sample_config["pipeline"] # pipeline/analysis to be executed
    # this stopped to --> workgetattr(__import__(command), "run")
    command_fn = getattr(globals()[pipeline], "run")
    command_fn(global_config, sample_config)


def check_consistency(global_config, sample_config):
    if "pipeline" not in sample_config:
        sys.exit("Error: pipeline must be specified in yaml sample "
                "configuration")
    if "genomeSize" not in sample_config:
        sys.exit("Error: genomeSize must be specified in yaml sample "
                "configuration (at least an estimate must be provided)")
    if "libraries" not in sample_config:
        sys.exit("Error: libraries must be specified in yaml sample "
                "configuration. At least one library must be specified")
    #check that the sepcified pipeline is supported
    pipelines = global_config["Pipelines"]
    pipeline  = sample_config["pipeline"]
    if pipeline not in pipelines:
        sys.exit("Error: pipeline {} does not exist. Supported pipelines "
                "are: {}".format(pipeline, pipelines))

    libraries=[library for library in sample_config["libraries"]]
    if len(libraries) < 1:
        sys.exit("Error: 0 libraries specified: at least one library must be "
                "present in the sample configuration file")
    for library, libraryData in sample_config["libraries"].items():
        if "pair1" not in libraryData or "pair2" not in libraryData \
            or "orientation" not in libraryData or "insert" not in \
            libraryData or "std" not in libraryData:
            sys.exit("Error: processing library {} each library must contatin "
                    "5 items: pair1, pair2, orientation, insert," 
                    "std".format(library))
        if not os.path.exists(libraryData["pair1"]):
            sys.exit("Error: pair1 {} does not exist".format(
                libraryData["pair1"]))
        if libraryData["pair2"] is not None:
            if not os.path.exists(libraryData["pair2"]):
                sys.exit("Error: pair2 {} does not exist".format(
                    libraryData["pair2"]))
        if libraryData["orientation"] not in ["innie", "outtie", "none"]:
            sys.exit("Error: orientation must be one among innie "
                    "(i.e., -> <-), outtie (i.e., <- ->), or none (single end)")
        if not isinstance(libraryData["insert"], int):
            sys.exit("Error: insert must be a integer")
        if not isinstance(libraryData["std"], int):
            sys.exit("Error: insert must be a integer")
#TODO: check that tools specified in the configuraiton are really working


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--global-config', help=("configuration file "
            "containing paths to tools"), type=str)
    parser.add_argument('--sample-config', help=("configuration file "
            "containing sample information"), type=str)
    args = parser.parse_args()

    main(args)

