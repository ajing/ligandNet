###############################
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

__AMINOACIDDIR__ = "/users/ajing/perllib/amino_acid_counter.pl"
#__BIOUNITDIR__ = "/users/ajing/BindingMOAD/biou/"
__BIOUNITDIR__ = "/users/ajing/BindingMOAD/biou_test/"
#__TMPDIR__ = "/users/ajing/ligandNet/pylib/BindingMoad2011/"
__EXTRACTDIR__ = "/users/ajing/ligandNet/pylib/BindingMoad2011_test/"
#__INPUTDIR__ = __EXTRACTDIR__ + "BindingMoad2011/"
__INPUTDIR__ = __EXTRACTDIR__ + "BindingMoad2011/"
#__WORKDIR__ = "/users/ajing/ligandNet/pylib/tmp/"
__WORKDIR__ = "/users/ajing/ligandNet/pylib/tmp_test/"

def callRichsCode(filedir, maxdistance, outputdir):
    commandlist = [ "perl", __AMINOACIDDIR__, filedir, "1", "notimportant.txt", maxdistance, outputdir]
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

def goThroughAllFiles():
    tmp_dir = __WORKDIR__ + "tmp.txt"
    final_dir = __WORKDIR__ + "final.txt"
    inputdir = __INPUTDIR__
    #try:
    #    shutil.rmtree(__TMPDIR__)
    #except:
    #    pass
    for filename in os.listdir( __BIOUNITDIR__ ):
        filedir = __BIOUNITDIR__ + filename
        unzipToTmp(filedir)
        inputfile = inputdir + os.listdir(inputdir)[0]
        for each in os.listdir(inputdir):
            if each.startswith(filename[:4]):
                inputfile = inputdir + each
                rich_out = __WORKDIR__ + each + "_out"
                callRichsCode( inputfile, str(4), rich_out)
                changeFileFormate( rich_out, tmp_dir )
                concatenateFiles( tmp_dir, final_dir )

if __name__ == "__main__":
    ## Test case
    #unzipToTmp("/users/ajing/BindingMOAD/1hg0.zip")
    #concatenateFiles( "tmp.txt", "final.txt")
    ## main function
    goThroughAllFiles()
