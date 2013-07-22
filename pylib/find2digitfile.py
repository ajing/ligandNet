"""
    This file is primarily for checking biounit files for whether 2 digit column name is very common and remove all BioUnit file with 2 digit chain id
    Author: ajing
    Data  : 7/17/2013
    Output: BioUnitWithtwoDigitChainID.txt,
"""
import os

BIOUNIT_DIR = "/users/ajing/ligandNet/pylib/BindingMoad2011_test/BindingMoad2011"
RESULT_DIR  = "/users/ajing/ligandNet/pylib/tmp_test"

def twoDigitChainID( biounit ):
    for line in open(biounit):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            ChainIDColumn = line[20:22]
            ChainIDColumnSplit = ChainIDColumn.split(" ")
            if len(ChainIDColumnSplit) == 1:
                return True
    return False

def returnAllChains( biounit ):
    allChains = []
    for line in open(biounit):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            ChainIDColumn = line[20:22]
            ChainIDColumnSplit = ChainIDColumn.split(" ")
            ChainID = ChainIDColumnSplit[-1]
            if ChainID == "z":
                continue
            if len(ChainID) > 1:
                ChainID = ChainID + "-"
            if not ChainID in allChains: # Here I assume only two situations: ["","A"] or ["AF"]
                allChains.append( ChainID )
    return allChains

def printAllBiounitHastwoDigitChainID():
    for filename in os.listdir( BIOUNIT_DIR ):
        filedir = os.path.join( BIOUNIT_DIR, filename)
        if twoDigitChainID( filedir ):
            os.remove( filedir )
        #ChainList = returnAllChains( filedir )
        #print filename + "\t" + ".".join(ChainList)

def notProcessYet( biounitfile ):
    for filename in os.listdir( RESULT_DIR ):
        if biounitfile[:4].lower() == filename[:4].lower():
            return False
    return True

def removeNotExecutiveFile():
    for filename in os.listdir( BIOUNIT_DIR ):
        filedir = os.path.join( BIOUNIT_DIR, filename)
        if notProcessYet( filename ):
            os.remove( filedir )

def debug():
    biounitdir = BIOUNIT_DIR + "/1a8s.bio2"
    print biounitdir

if __name__ == "__main__":
    #debug()
    #removeNotExecutiveFile()
    printAllBiounitHastwoDigitChainID()
