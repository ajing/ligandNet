"""
    This program is majorly for get the chain list for protein, and check whether a ligand is in bindingMOAD list or not.
    No chain z and no chain with length less than 10
    Author: ajing
    Date:   7/17/2013
"""

## For every_parser
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

from Modular import BIOUNIT_DIR, EVERYCSV

import Bio
from Bio.PDB.PDBExceptions import PDBConstructionException
############## exception for ligand which is less than 10 residues, but not in every.csv file ##############
__EXCEPTION_FILE__ = "/users/ajing/ligandNet/Error/proteinChainException.txt"
__ERROR_FILE__ = "/users/ajing/ligandNet/Error/proteinChainError.txt"
__BIOUNIT_DIR__ = BIOUNIT_DIR
__BIOUNIT_CANNOT_READ__ = "defaultProteinChainList.txt"
__OUTPUT_PROTEINLIGAND__ = "proteinChain.txt"

from multiprocessing import Pool
# a global variable for ALL
__ALL__ = dict()

# For getting files in directory:
import os

def getChainLength( chain ):
    residuenum = 0
    for residue in chain.get_list():
        residuenum = residuenum + 1
    return residuenum

def makeStructure(filename):
    """
        return structure for a file
    """
    #pdb_code = filename.split(".")[0]
    pdb_code = filename.split("/")[-1]
    pdb_filename = filename
    error_obj = open( __ERROR_FILE__, "a" )
    try:
        structure = Bio.PDB.PDBParser(PERMISSIVE = 1).get_structure(pdb_code, pdb_filename)
    except PDBConstructionException:
        error_obj.write( "some value error when reading: " + pdb_code + "\n")
        error_obj.close()
        return False
    except:
        error_obj.write( "some value error when reading: " + pdb_code + "\n")
        error_obj.close()
        return False
    error_obj.close()
    return structure

def residueNameinLigandlist( residuename, ligandlist ):
    # Deal with crazy long peptide name
    if residuename in ligandlist:
        return True
    for each in ligandlist:
        if residuename in each:
            return True
    return False

def returnProteinChainID( structure, filename, ligandlist ):
    proteinChains = []
    exception_obj = open( __EXCEPTION_FILE__, "a")
    for model in structure.get_list():
        for chain in model.get_list():
            if chain.id == "z":
                continue
            if getChainLength( chain ) > 10 and not chain.id in proteinChains:
                proteinChains.append( chain.id )
            else:
                for residue in chain.get_list():
                    residueName = residue.id[0]
                    if residueName.startswith("H_"):
                        residueName = residueName[2:]
                    if not residueNameinLigandlist( residueName, ligandlist ):
                        exception_obj.write( "Cannot Find Ligand in MOAD:" + "\t".join( [filename, residueName, chain.id] ) + "\n" )
    return proteinChains

def returnProteinChainIDSimple( biounitdir ):
    # this function only take SEQRES for all chain ids
    # For simplicity here we only keep chain id with more than 13 residues.
    exception_obj = open( __EXCEPTION_FILE__, "a")
    exception_obj.write("Cannot Find Protein Chain For:" + biounitdir)
    chainIDList = []
    proteinChains = []
    for line in open( biounitdir ):
        if line.startswith( "SEQRES" ):
            chainID = line[11]
            if chainID in chainIDList and not chainID in proteinChains:
                proteinChains.append( chainID )
            else:
                chainIDList.append( chainID )
    proteinChains.sort()
    return proteinChains

def processOneBioUnit( biounitdir ):
    out_obj = open( __OUTPUT_PROTEINLIGAND__, "a")
    filename  = biounitdir.split("/")[-1]
    PDBID     = filename[:4].upper()
    try:
        structure = makeStructure( biounitdir )
        PDBChainIDs = returnProteinChainID( structure, filename, __ALL__[PDBID].keys() )
    except:
        print "something wrong with " + filename
        PDBChainIDs = returnProteinChainIDSimple( biounitdir )
    out_obj.write( filename + "\t" + "".join(PDBChainIDs) + "\n" )
    out_obj.close()

def BioUnitFilter( BioUnitDirectory, existingfile ):
    # return new biounit file list without the file already in __OUTPUT_PROTEINLIGAND__
    filedir = os.listdir( BioUnitDirectory )
    # existing protein id don't need to process
    proteinIDList = []
    PDBList = []
    for line in open(existingfile):
        content = line.strip().split("\t")
        proteinIDList.append( content[0].lower() )
        PDBID   = content[0].split(".")[0]
        if not PDBID in PDBList:
            PDBList.append( PDBID )
    NewDirectoryList = []
    totalPDBBioUnitDirectory = []
    for each in filedir:
        PDB = each.split(".")[0]
        if not PDB in totalPDBBioUnitDirectory:
            totalPDBBioUnitDirectory.append( PDB )
        if not each.lower() in proteinIDList:
            NewDirectoryList.append( each )
    print "total number of files in : " + BioUnitDirectory
    print len(filedir)
    print "Distinctive PDBs in : " + BioUnitDirectory
    print len(totalPDBBioUnitDirectory)
    print "Distinctive PDBs in existing file : " + existingfile
    print len(PDBList)
    print "BioUnit files I still work on for : " + BioUnitDirectory
    print len(NewDirectoryList)
    return NewDirectoryList

def runOneBioUnit( biounitdirlist ):
    error_obj = open( __ERROR_FILE__, "a" )
    for each in biounitdirlist:
        try:
            processOneBioUnit( each )
        except:
            print "something wrong with " + each
            error_obj.write( "something wrong with " + each + "\n" )
    error_obj.close()

def callbackforBioUnit( biounitdir ):
    error_obj = open( __ERROR_FILE__, "a" )
    error_obj.write( "something wrong with " + biounitdir + "\n" )
    error_obj.close()

def removeExceptionFile( BioUnitDirList ):
    fileCannotRead = []
    for line in open(__BIOUNIT_CANNOT_READ__):
        content = line.strip().split("\t")
        fileCannotRead.append( content[0] )
    newList = []
    for each in BioUnitDirList:
        if not each in fileCannotRead:
            newList.append( each )
    return newList

def addDefaultBioUnit():
    # some biounit files cannot be processed by biopython, so here I just add it as default using information in "SEQRES"
    output_obj = open( __OUTPUT_PROTEINLIGAND__, "a")
    for line in open(__BIOUNIT_CANNOT_READ__):
        output_obj.write( line )
    output_obj.close()

def main():
    # remove error and output files first
    try:
        #os.remove( __ERROR_FILE__ )
        os.remove( __EXCEPTION_FILE__ )
    except:
        pass
    if not os.path.exists(__OUTPUT_PROTEINLIGAND__):
        addDefaultBioUnit()
    everyparser = every_parser(EVERYCSV)
    everyparser.find_PDBID_ValidLigand()
    global __ALL__
    __ALL__ = everyparser.ALL
    pool = Pool(processes = 7)
    argumentlist = [ __BIOUNIT_DIR__ + filename for filename in removeExceptionFile( BioUnitFilter(__BIOUNIT_DIR__, __OUTPUT_PROTEINLIGAND__) ) ]
    ####################### 8/19/2013 checking for what's wrong with each file  #########################
    runOneBioUnit(argumentlist)
    #####################################################################################################
    print argumentlist
    #result = pool.map_async( processOneBioUnit, argumentlist )
    #resulttxt = result.wait()
    print resulttxt

if __name__ == "__main__":
    main()
