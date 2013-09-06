##############################
##   Date  : 06/18/2013
##   Author: ajing
##   Purpose: parse ligand interaction output from Rich
##            And prepare database input for BindingPocketInfo
##   I/O   : input is the output of amino_acid_counter.pl
##           And output is a concatenation file of all output of amino_acid_counter.pl
###############################

import zipfile
import os
import subprocess
import shutil
from multiprocessing import Pool
from functools import partial
import time
import random
# a timeout wrapper for function
import timeout

__AMINOACIDDIR__ = "/users/dicksmit/amino_acid_counter.pl"
__BIOUNITDIR__ = "/users/ajing/BindingMOAD/biou_test/"
__INPUTDIR__ = "/users/ajing/ligandNet/2012_biounits/"
__WORKDIR__ = "/users/ajing/ligandNet/tmpfiles/"

#################### Parameter ########################
__MAX_DISTANCE__ = 4.5

def callRichsCode(filedir, outputdir):
    maxdistance = str(__MAX_DISTANCE__)
    commandlist = [ "perl", __AMINOACIDDIR__, filedir, "2", "notimportant2.txt", maxdistance, outputdir]
    print commandlist
    command = subprocess.list2cmdline( commandlist)
    cm = subprocess.Popen( command, shell = True, stdout = subprocess.PIPE ).communicate()

def concatenateFiles(orifile, destfile):
    destination = open( destfile, "a")
    original = open(orifile, "rb")
    destination.write( original.read() )
    destination.close()

def changeFileFormate( filedir, tmpdir ):
    # add first column as PDBID and strip some space
    PDBID = filedir.split("/")[-1].split("_")[-2]
    print PDBID
    tmpobj = open( tmpdir, "w")
    for line in open(filedir):
        content = line.strip().split("\t")
        newcontent = []
        for each in content:
            newcontent.append(each.strip())
        newcontent.pop(5)
        newcontent.insert( 0, PDBID.upper() )
        if newcontent[-1] == "Metal":
            newcontent[-1] = str(1)
        else:
            newcontent.append( str(0) )
        tmpobj.write( ",".join(newcontent) + "\n")
    tmpobj.close()

def partial_operateOneFileWithExistingFiles( filedir ):
    # Here is the modified code which just use the files under __INPUTDIR__
    final = __WORKDIR__ + "final.txt"
    inputdir = __INPUTDIR__
    workdir = __WORKDIR__
    filename = filedir.split("/")[-1]
    inputfile = inputdir + os.listdir(inputdir)[0]
    rich_out = workdir + filename + "_out"
    tmp = workdir + filename + "_tmp"
    callRichsCode( filedir, rich_out )    # str(4) means the maximun distance is 4, also mentioned in callRichsCode
    changeFileFormate( rich_out, tmp )
    concatenateFiles( tmp, final )

def filterInputFile(filename, inputdir):
    inputlist = os.listdir(inputdir)
    inputlistrefine = [ i.split(".")[0] for i in inputlist]
    if filename.split(".")[0] in inputlistrefine:
        return False
    return True

def goThroughAllFiles():
    # For only run on file in __INPUTDIR__
    pool = Pool(processes = 6)
    inputdir = __INPUTDIR__
    argumentlist = [ __INPUTDIR__ + filename for filename in os.listdir( __INPUTDIR__ ) ]
   # print argumentlist
    result = pool.map_async( partial_operateOneFileWithExistingFiles, argumentlist )
    resulttxt = result.get()
    resultfile = open("poolresult.txt", "w")
    resultfile.write(resulttxt)
    #    print result.get()

if __name__ == "__main__":
    goThroughAllFiles()
