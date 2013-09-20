import sys, os, yaml, glob
import argparse
import QCcontrol


def main(args):
    with open(args.global_config) as in_handle:
        global_config = yaml.load(in_handle)

    with open(args.sample_config) as sample_config_handle:
        sample_config = yaml.load(sample_config_handle)

    check_consistency(global_config,sample_config)

    #TODO: prepare the folder structure with softlins and change the location of the fails to point to the softlink position
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    if os.path.exists(os.path.join(os.getcwd(), "DATA")):
        sorted_libraries_by_insert = update_sample_config(sorted_libraries_by_insert)
    else:
        sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)
    
    run_analys(global_config, sample_config, sorted_libraries_by_insert)
    

    """for operation, instruments in global_config.items():
        print "{}".format(operation)
        for tool, tool_info in instruments.items():
            print "\t{}".format(tool)
            print "\t\t{}".format(tool_info["bin"])
    """

def prepare_folder_structure(sorted_libraries_by_insert):
    mainDir = os.getcwd()
    DataFolder = os.path.join(os.getcwd(), "DATA")
    if os.path.exists(DataFolder):
        sys.exit("DATA dir already exists: danger to over-write data: terminate execution")
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
        pair1, pair2 = createSoftLinks(pair1, pair2, orientation, type, currentLibraryNumber)
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
        libraryInfo["pair1"] = _new_name(pair1, orientation, type, currentLibraryNumber, 1)
        libraryInfo["pair2"] = _new_name(pair2, orientation, type, currentLibraryNumber, 2)
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

def _new_name(oldPathName, orientation, type, currentLibraryNumber, pairNumber):
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


def run_analys(global_config, sample_config, sorted_libraries_by_insert):
    analysis = sample_config["operation"].keys()[0] # analysis to be executed
    command = analysis
    command_fn = getattr(__import__(command), "run")
    command_fn(global_config, sample_config, sorted_libraries_by_insert)
    


#def analysis_initiate(global_config, sample_config):
#    print "going to evalaute input data for de novo"


def check_consistency(global_config, sample_config):
    if "operation" not in sample_config:
        sys.exit("Error: operations must be specified in yaml sample configuration")
    if "genomeSize" not in sample_config:
        sys.exit("Error: genomeSize must be specified in yaml sample configuration (at least an estimate must be provided)")
    if "kmer" not in sample_config:
        sys.exit("Error: kmer must be specified in yaml sample configuration")
    if "libraries" not in sample_config:
        sys.exit("Error: libraries must be specified in yaml sample configuration. At least one library must be specified")
    #check that a possible opration is specified
    operations = [operation for operation in global_config]
    operation = sample_config["operation"].keys()[0]
    if operation not in operations:
        sys.exit("Error: operation {} does not exist. Supported operations are: {}".format(operation, operations))
    tool = sample_config["operation"][operation]["tool"]
    supportedTools = global_config[operation].keys()
    if tool not in supportedTools:
        sys.exit("Error: tool {} not supported for operation {}. Only the following tools can be used: {}".format(tool, operation, supportedTools))
    libraries=[library for library in sample_config["libraries"]]
    if len(libraries) < 1:
        sys.exit("Error: 0 libraries specified: at least one library must be present in the sample configuration file")
    for library, libraryData in sample_config["libraries"].items():
        if "pair1" not in libraryData or "pair2" not in libraryData \
            or "orientation" not in libraryData or "insert" not in libraryData or "std" not in libraryData:
            sys.exit("Error: processing library {} each library must contatin 5 items: pair1, pair2, orientation, insert, std".format(library))
        if not os.path.exists(libraryData["pair1"]):
            sys.exit("Error: pair1 {} does not exist".format(libraryData["pair1"]))
        if libraryData["pair2"] is not None:
            if not os.path.exists(libraryData["pair2"]):
                sys.exit("Error: pair2 {} does not exist".format(libraryData["pair2"]))
        if libraryData["orientation"] not in ["innie", "outtie", "none"]:
            sys.exit("Error: orientation must be one among innie (i.e., -> <-), outtie (i.e., <- ->), or none (single end)")
        if not isinstance(libraryData["insert"], int):
            sys.exit("Error: insert must be a integer")
        if not isinstance(libraryData["std"], int):
            sys.exit("Error: insert must be a integer")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--global-config', help="configuration file containing paths to tools", type=str)
    parser.add_argument('--sample-config', help="configuration file containing sample information", type=str)
    args = parser.parse_args()

    main(args)