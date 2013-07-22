"""
    Purpose: checking the consistent between BindingMOAD and Biounit files
    Date:    07/01/2013
    Author:  ajing
"""

## For every_parser
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionException
import os
__BIOUNIT_UNZIP_DIR__ = "/users/ajing/BindingMOAD/biou_unzip/"
## For testing
#__BIOUNIT_UNZIP_DIR__ = "/users/ajing/BindingMOAD/biou_part/"
## warning catch system is not really working, so __ERROR_FILE__ is not really useful here
__ERROR_FILE__ = "/users/ajing/ligandNet/pylib/compareError.txt"
__COMPARE_FILE__ = "/users/ajing/ligandNet/pylib/compare.txt"
__VALID_LIGAND_FILE__ = "/users/ajing/ligandNet/pylib/ligand_valid.txt"
__INVALID_LIGAND_FILE__ = "/users/ajing/ligandNet/pylib/ligand_invalid.txt"
__EXCEPTION_FILE__ = "/users/ajing/ligandNet/pylib/compareException.txt"

from multiprocessing import Pool
# a global variable for ALL
__ALL__ = dict()

from timeout import *

def printAllPossibleLowercaseChain( filedir ):
    """
       I want to check whether people usually use lowercase to mark ligands, but actually chain name for ligands can be anything
    """
    filelist = os.listdir( filedir )
    filedirlist = [ filedir + i for i in filelist]
    for file in filedirlist:
        structure = makeStructure( file )
        lowerCase = isLowerCaseChain( structure )
        if lowerCase:
            print lowerCase

def makeStructure(filename):
    """
        return structure for a file
    """
    #pdb_code = filename.split(".")[0]
    pdb_code = filename.split("/")[-1]
    pdb_filename = filename
    error_obj = open( __ERROR_FILE__, "w" )
    try:
        structure = Bio.PDB.PDBParser(PERMISSIVE = 1).get_structure(pdb_code, pdb_filename)
    except PDBConstructionException:
        error_obj.write( w.message )
    except:
        error_obj.write( "some value error when reading: " + pdb_code )
        return False
    error_obj.close()
    return structure

def isLowerCaseChain( structure ):
    """
        return chain id if chain id is lowercase
    """
    for model in structure.get_list():
        for chain in model.get_list():
            if chain.id != 'z' and chain.id.islower():
                return chain.id
    return False

def ligandListFromStruct( structure ):
    """
        structure is Bio.PDB structure, return a set of 3 letter string for ligands
        Here we will return a ligand list without "H_"
    """
    ligandlist = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                #print residue.id
                residueName = residue.id[0]
                if residueName.startswith("H_"):
                    ligandlist.append(residueName[2:].strip())
    return list(set(ligandlist))

def ligandListBindingMOAD( PDBID, ALL ):
    return [ i.strip() for i in ALL[PDBID].keys() ]

def diffTwoLists( list1, list2 ):
    set1  = set( list1 )
    set2  = set( list2 )
    list1 = list( set1.difference(set2) )
    list2 = list( set2.difference(set1) )
    overlap = list( set1 & set2 )
    return ( list1, list2, overlap )

def printHETResidue( residue ):
    """
        return line format as Aqeel told me.
    """
    residueID   = residue.get_full_id()
    proteinID   = residueID[0]
    chainID     = residueID[2]
    residueName = residueID[3][0]
    residueNum  = residueID[3][1]
    resultlist  = [ proteinID, residueName, residueNum, chainID]
    return "\t".join( [ str(i) for i in resultlist ] )

def printAllResiduewithLigandName( ligandList, structure ):
    """
        from a list of hetero residue list, return the output string for final writing to file
    """
    outputString = ""
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                residueName = residue.id[0]
                if residueName.startswith("H_") and residueName[2:] in ligandList:
                    outputString = outputString + printHETResidue(residue) + "\n"
    return outputString

def returnValidInvalid( PDBID, ligandList, ALL ):
    validlist = []
    invalidlist = []
    for ligand in ligandList:
        if ALL[PDBID][ligand]:
            validlist.append( ligand )
        else:
            invalidlist.append( ligand )
    return (validlist, invalidlist)

def returnProperBIOUNIT(PDBID, ligandList ):
    """
        Here I assume there is one biounit contains all valid ligands
    """
    filelist = os.listdir( __BIOUNIT_UNZIP_DIR__ )
    filedirlist = [ __BIOUNIT_UNZIP_DIR__ + filename for filename in filelist if filename.startswith( PDBID.lower() )]
    for eachfile in filedirlist:
        structure   = makeStructure(eachfile)
        if not structure:
            continue
        BiounitList = ligandListFromStruct( structure )
        if set(ligandList).issubset( set(BiounitList) ) or set(ligandList) == set(BiounitList):
            return eachfile
    return False

def returnProperBIOUNIToneLigand(PDBID, aligand ):
    """
        Here I assume there is one biounit contains all valid ligands
    """
    filelist = os.listdir( __BIOUNIT_UNZIP_DIR__ )
    filedirlist = [ __BIOUNIT_UNZIP_DIR__ + filename for filename in filelist if filename.startswith( PDBID.lower() )]
    for eachfile in filedirlist:
        structure   = makeStructure(eachfile)
        if not structure:
            continue
        BiounitList = ligandListFromStruct( structure )
        if aligand in BiounitList:
            return eachfile
    return False

@timeout(10)
def processOnePDBID( eachPDBID ):
    valid_ligand_file = open( __VALID_LIGAND_FILE__, "a" )
    noSuchBiounit = open( __EXCEPTION_FILE__, "a")
    for eachligand in __ALL__[ eachPDBID ].keys():
        Biounitfile = returnProperBIOUNIToneLigand( eachPDBID, eachligand)
        if Biounitfile:
            structure   = makeStructure(Biounitfile)
            valid_ligand_file.write( printAllResiduewithLigandName( [eachligand], structure ) )
        else:
            noSuchBiounit.write( "cannot find proper biounit for " + eachPDBID + "\n")
    valid_ligand_file.close()
    noSuchBiounit.close()

def mainCompare2():
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    global __ALL__
    __ALL__ = everyparser.ALL
    pool = Pool(processes = 6)
    argumentlist = __ALL__.keys()
    print argumentlist
    result = pool.map_async( processOnePDBID, argumentlist)
    try:
        resulttxt = result.get()
    except TimeoutError as e:
        print e.error_message

    #print dir(resulttxt)
    #resultfile = open("poolresult.txt", "w")
    #resultfile.write(resulttxt)
    #for eachPDBID in __ALL__.keys():
    #    processOnePDBID( eachPDBID )

def mainCompare():
    """
        main function for comparing biounit and MOAD list. Assumption here is that system should be linux with "/", but  all the directories are under linux, so should be fine.
    """
    compare_file = open( __COMPARE_FILE__, "w" )
    valid_ligand_file = open( __VALID_LIGAND_FILE__, "w" )
    invalid_ligand_file = open( __INVALID_LIGAND_FILE__, "w" )
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL = everyparser.ALL
    filelist = os.listdir( __BIOUNIT_UNZIP_DIR__ )
    filedirlist = [ __BIOUNIT_UNZIP_DIR__ + filename for filename in filelist ]
    #filedirlist = [ __BIOUNIT_UNZIP_DIR__ + "1in5.bio1" ]
    for eachfile in filedirlist:
        structure   = makeStructure( eachfile )
        # pass the exception of error biounit file
        if not structure:
            continue
        BiounitList = ligandListFromStruct( structure )
        print "BiounitList"
        print BiounitList
        PDBID       = eachfile.split("/")[-1].split(".")[0].upper()
        MOADList    = ligandListBindingMOAD( PDBID, ALL )
        print "MOADList"
        print MOADList
        BiounitMore, MOADMore, Overlap = diffTwoLists( BiounitList, MOADList )
        print (BiounitMore, MOADMore, Overlap)
        eachfilename = eachfile.split("/")[-1]
        if BiounitMore:
            compare_file.write( stringDiff( BiounitMore, "BIOUNIT " + eachfilename + " compare to " + PDBID ) + "\n" )
        if MOADMore:
            compare_file.write( stringDiff( MOADMore, "MOAD " + PDBID + " compare to " + eachfilename) + "\n" )
        if Overlap:
            validlist, invalidlist = returnValidInvalid( PDBID, Overlap, ALL)
            valid_ligand_file.write( printAllResiduewithLigandName( validlist, structure ) )
            invalid_ligand_file.write( printAllResiduewithLigandName( invalidlist, structure ) )
    compare_file.close()
    valid_ligand_file.close()
    invalid_ligand_file.close()


def stringDiff( MoreList, type):
    """
        Return string for list from diffTwoLists
    """
    return type + " has more ligands: " + ",".join( MoreList )

def debug():
    ## Debug for everyparser
    #everyparser = every_parser()
    #everyparser.find_PDBID_ValidLigand()
    #print everyparser.ALL
    ## Debug for PDB reader
    os.chdir(__BIOUNIT_UNZIP_DIR__)
    pdb_code = "1XI4"
    pdb_filename = "1glg.bio1"
    structure = Bio.PDB.PDBParser(PERMISSIVE = 1).get_structure(pdb_code, pdb_filename)
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if chain.id == 'z':
                    residueid = residue.get_full_id()
                    print residue.get_full_id()
                    print residueid[3][0]
            #print dir(chain)
            #print chain.get_id()
            #print chain.get_full_id()
            #print chain.get_level()
            #print chain.xtra


if __name__ == "__main__":
    #debug()
    #printAllPossibleLowercaseChain( __BIOUNIT_UNZIP_DIR__ )
    # Test
    #print diffTwoLists( [2,1,3,4], [2,3])
    os.remove( __VALID_LIGAND_FILE__ )
    os.remove( __EXCEPTION_FILE__ )
    mainCompare2()
