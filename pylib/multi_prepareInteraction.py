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
#from timeout import timeout

#__AMINOACIDDIR__ = "/users/ajing/perllib/amino_acid_counter.pl"
__AMINOACIDDIR__ = "/users/dicksmit/amino_acid_counter.pl"
#__BIOUNITDIR__ = "/users/ajing/BindingMOAD/biou/"
__BIOUNITDIR__ = "/users/ajing/BindingMOAD/biou_test/"
#__TMPDIR__ = "/users/ajing/ligandNet/pylib/BindingMoad2011/"
__EXTRACTDIR__ = "/users/ajing/ligandNet/pylib/BindingMoad2011_test/"
#__INPUTDIR__ = __EXTRACTDIR__ + "BindingMoad2011/"
__INPUTDIR__ = __EXTRACTDIR__ + "BindingMoad2011/"
#__WORKDIR__ = "/users/ajing/ligandNet/pylib/tmp/"
__WORKDIR__ = "/users/ajing/ligandNet/pylib/tmp_test/"

#################### Parameter ########################
__MAX_DISTANCE__ = 4.5

#@timeout(10)
def callRichsCode(filedir, outputdir):
    maxdistance = str(__MAX_DISTANCE__)
    commandlist = [ "perl", __AMINOACIDDIR__, filedir, "2", "notimportant2.txt", maxdistance, outputdir]
    print commandlist
    command = subprocess.list2cmdline( commandlist)
    cm = subprocess.Popen( command, shell = True, stdout = subprocess.PIPE ).communicate()

def unzipToTmp(filedir):
    myzip = zipfile.ZipFile(filedir, "r")
    myzip.extractall(__EXTRACTDIR__)

def concatenateFiles(orifile, destfile):
    destination = open( destfile, "a")
    original = open(orifile, "rb")
    destination.write( original.read() )
    destination.close()
    # remove the orifile
    #shutil.rmtree(__TMPDIR__)

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

def operateOneFile( filedir, inputdir, workdir, final ):
    filename = filedir.split("/")[-1]
    unzipToTmp(filedir)
    inputfile = inputdir + os.listdir(inputdir)[0]
    for each in os.listdir(inputdir):
        if each.startswith(filename[:4]):
            inputfile = inputdir + each
            rich_out = workdir + each + "_out"
            tmp = workdir + each + "_tmp"
            callRichsCode( inputfile, rich_out)    # str(4) means the maximun distance is 4, also mentioned in callRichsCode
            changeFileFormate( rich_out, tmp )
            concatenateFiles( tmp, final )

def partial_operateOneFile( filedir ):
    # 6/26 Because there is not good functool partial function, so I wrote this ugly function for partially argument of operateOneFile
    final = __WORKDIR__ + "final.txt"
    inputdir = __INPUTDIR__
    workdir = __WORKDIR__
    filename = filedir.split("/")[-1]
    unzipToTmp(filedir)
    inputfile = inputdir + os.listdir(inputdir)[0]
    for each in os.listdir(inputdir):
        if each.startswith(filename[:4]):
            inputfile = inputdir + each
            rich_out = workdir + each + "_out"
            tmp = workdir + each + "_tmp"
            callRichsCode( inputfile, rich_out)    # str(4) means the maximun distance is 4, also mentioned in callRichsCode
            changeFileFormate( rich_out, tmp )
            concatenateFiles( tmp, final )

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
    pool = Pool(processes = 8)
    inputdir = __INPUTDIR__
    argumentlist = [ __BIOUNITDIR__ + filename for filename in os.listdir( __BIOUNITDIR__ ) if filterInputFile(filename, inputdir) ]
    print argumentlist
    result = pool.map_async( partial_operateOneFile, argumentlist )
    resulttxt = result.get()
    #for filename in os.listdir( __BIOUNITDIR__ ):
    #    filedir = __BIOUNITDIR__ + filename
    #    #p = multiprocessing.Process( target = operateOneFile, args = (filedir, inputdir, __WORKDIR__, final_dir))
    #    #p.start()
    #    result = pool.apply_async(operateOneFile, (filedir, inputdir, __WORKDIR__, final_dir))

def goThroughAllFiles2():
    # For only run on file in __INPUTDIR__
    pool = Pool(processes = 8)
    inputdir = __INPUTDIR__
    argumentlist = [ __INPUTDIR__ + filename for filename in os.listdir( __INPUTDIR__ ) ]
   # print argumentlist
    result = pool.map_async( partial_operateOneFileWithExistingFiles, argumentlist )
    resulttxt = result.get()
    resultfile = open("poolresult.txt", "w")
    resultfile.write(resulttxt)
    #    print result.get()

if __name__ == "__main__":
    ## Test case
    #unzipToTmp("/users/ajing/BindingMOAD/1hg0.zip")
    #concatenateFiles( "tmp.txt", "final.txt")
    ## main function
    #goThroughAllFiles()
    goThroughAllFiles2()
