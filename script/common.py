import sys, os, yaml, glob
import subprocess

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



def directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        return 0
    else:
        print "done ({} folder already present, assumed already run)".format(directory)
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
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    if os.path.exists(os.path.join(os.getcwd(), "DATA")):
        sorted_libraries_by_insert = update_sample_config(sorted_libraries_by_insert)
    else:
        sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)
    return sorted_libraries_by_insert



