from __future__ import absolute_import
from __future__ import print_function
import sys, os, yaml, glob
import subprocess
import datetime
import time
import types


# Header levels
H1, H2, H3, H4, H5, H6 = 1, 2, 3, 4, 5, 6 

# List styles
UL, OL = 0, 1

# Alignment
CENTER, LEFT, RIGHT = 'CENTER', 'LEFT', 'RIGHT'


def prepare_folder_structure(sorted_libraries_by_insert):
    mainDir = os.getcwd()
    DataFolder = os.path.join(os.getcwd(), "DATA")
    if os.path.exists(DataFolder):
        sys.exit("DATA dir already exists: danger to over-write data: \
                terminate execution")
    os.makedirs(DataFolder)
    os.chdir(DataFolder)
    CurrentDir = os.getcwd()
    #now prepare softlinks to data and give to libraries human readable names
    currentLibraryNumber = 1;
    type = ["SE", "PE", "MP"]
    for library, libraryInfo in sorted_libraries_by_insert:
        pair1 = libraryInfo["pair1"]
        pair2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        pair1, pair2 = createSoftLinks(pair1, pair2, orientation, type,
                currentLibraryNumber)
        libraryInfo["pair1"] = pair1
        libraryInfo["pair2"] = pair2
        currentLibraryNumber += 1
    os.chdir("..")
    return sorted_libraries_by_insert

def update_sample_config(sorted_libraries_by_insert):
    mainDir = os.getcwd()
    DataFolder = os.path.join(os.getcwd(), "DATA")
    if not os.path.exists(DataFolder):
        sys.exit("DATA dir does not exists: we should not be here!!!!")
    os.chdir(DataFolder)
    currentLibraryNumber = 1
    type = ["SE", "PE", "MP"]
    for library, libraryInfo in sorted_libraries_by_insert:
        pair1 = libraryInfo["pair1"]
        pair2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        libraryInfo["pair1"] = _new_name(pair1, orientation, type,
                currentLibraryNumber, 1)
        libraryInfo["pair2"] = _new_name(pair2, orientation, type,
                currentLibraryNumber, 2)
        currentLibraryNumber += 1
    os.chdir("..")
    return sorted_libraries_by_insert

def createSoftLinks(pair1, pair2, orientation, type, currentLibraryNumber):
    pair1NewName = _new_name(pair1, orientation, type, currentLibraryNumber, 1)
    pair2NewName = _new_name(pair2, orientation, type, currentLibraryNumber, 2)
    os.symlink(pair1, pair1NewName)
    if pair2NewName is not None:
         os.symlink(pair2, pair2NewName)
    return pair1NewName, pair2NewName

def _new_name(oldPathName, orientation, type, currentLibraryNumber,
        pairNumber):
    if oldPathName is None:
        return oldPathName;
    oldName = os.path.split(oldPathName)[1]
    oldNameHead , oldNameTail = oldName.split(".",1)
    newName = "lib{}_".format(currentLibraryNumber)
    if orientation == "none":
        newName += "SE."
    elif orientation == "innie":
        newName += "PE_{}.".format(pairNumber)
    elif orientation == "outtie":
        newName += "MP_{}.".format(pairNumber)
    newName += oldNameTail
    newName = os.path.join(os.getcwd(), newName)
    return newName



def directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        return 0
    else:
        print("done ({} folder already present, assumed already run)".format(
                directory))
        return 1

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None



def _sort_libraries_by_insert(sample_config):
    sorted_libraries_by_insert = sorted(sample_config["libraries"].items(),
            key=lambda v: v[1]["insert"])
    return sorted_libraries_by_insert



def check_dryrun(sample_config):
    if "dryrun" in sample_config:
        return 1
    return 0

def get_command_str(command):
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    if type(command) is list:
        command = ' '.join(command)

    return "{}\n{}\n".format(st, command)

def print_command(command):
    """ prints in a human readable way a command stored in a list"""
    print(get_command_str(command))


def _check_pipeline(sample_config, global_config):
    """ check that user specified commands are supported by this pipeline"""
    print("checking tool consistency:")

    pipeline     = sample_config["pipeline"]
    user_tools   = sample_config["tools"] #might be empty
    global_tools = global_config["Pipelines"][pipeline]
    """be sure that pipeline I am going to run tries to execute only supported
    programs"""
    for tool in user_tools:
        if tool not in global_tools:
            print("tool {} not available in pipeline {}. Only these methods "
                   "are available: {}".format(tool, pipeline, global_tools))
            sys.exit("Error: invalid configuration file(s)")

    """check that programs to be executed are supported"""
    tools_to_check = []
    # check that the user wants to run only some (or all) the programs
    if user_tools is not []:
        tools_to_check = [i for i in user_tools]
    else:
        # otherwise check all the tools supported
        tools_to_check = [i for i in global_config["Pipelines"][pipeline]]
    if "align" in tools_to_check:
        tools_to_check.remove("align")
        tools_to_check.extend(["bwa", "samtools"])
        picard_dir = ""
        if os.environ.get('PICARD_HOME'):
            picard_dir = os.environ.get('PICARD_HOME')
        elif "picard" in global_config["Tools"]:
            picard_dir = global_config["Tools"]["picard"]["bin"]
        if not os.path.isdir(picard_dir):
            sys.exit("align is part of the pipeline you want to run. "
                    "PicardTools must be present in the enviorment or "
                    "specified in global config "
                    "directory {} does not exists.".format(picard_dir))


    for tool in tools_to_check:
        if tool not in global_config["Tools"]:
            sys.exit("Tool {} is not present as a valid entry in global "
                    "configuration file. Please add it in order to execute the "
                    "command.".format(tool))

        tool_bin = global_config["Tools"][tool]["bin"]
        print(tool_bin)
        """this step is a case by case step as several special cases are 
        present"""
        special_tools = ["allpaths", "abyss", "cabog", "masurca",  "trinity",
                "trimmomatic"]
        if tool in special_tools:
            binaries_to_check = []
            if tool == "allpaths":
                binaries_to_check.append(os.path.join(tool_bin,
                    "PrepareAllPathsInputs.pl"))
                binaries_to_check.append(os.path.join(tool_bin,
                    "RunAllPathsLG"))
            if tool == "abyss":
                binaries_to_check.append(os.path.join(tool_bin,
                    "abyss-pe"))
            if tool == "cabog":
                binaries_to_check.append(os.path.join(tool_bin,"runCA"))
            if tool == "masurca":
                binaries_to_check.append(os.path.join(tool_bin,"bin",
                    "runSRCA.pl"))
            if tool == "picard":
                binaries_to_check.append(os.path.join(tool_bin,))
            if tool == "trinity": ## check special tools
                binaries_to_check.append(os.path.join(tool_bin,
                    "Trinity"))
                binaries_to_check.append(os.path.join(tool_bin, "util",
                    "align_and_estimate_abundance.pl"))
            if tool == "trimmomatic":
                if not os.path.exists(tool_bin):
                    sys.exit("tool trimmomatic specified but file {} does not "
                            "exists. If you are not planning to run this tool "
                            "specify the list of to be run tools in the "
                            "sample_config section.".format(tool_bin))

            for binary in binaries_to_check:
                if not which(binary):
                    sys.exit("tool {} requires availability of binary {} but "
                            "pipeline did not find it. Please check the tool "
                            "installation. If you are not planning to run this "
                            "tool specify the list of to be run tools in the "
                            "sample_config section.".format(tool, binary))
        else:
            if not which(tool_bin):
                 sys.exit("Error: tool {} requires full path to binaries in "
                         "config file. Path {} does not exists. Please modify "
                         "the global config.If you are not planning to run "
                         "this tool specify the list of to be run tools in the "
                         "sample_config section.".format(tool, tool_bin))
    print("all tools the pipeline is going to use appear properly installed.")
    return

