"""
    This program is majorly for get the chain list for protein, and check whether a ligand is in bindingMOAD list or not.
"""

## For every_parser
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

import Bio
from Bio.PDB.PDBExceptions import PDBConstructionException
__EXCEPTION_FILE__ = "/users/ajing/ligandNet/pylib/proteinChainException.txt"
__ERROR_FILE__ = "/users/ajing/ligandNet/pylib/proteinChainError.txt"
__BIOUNIT_DIR__ = "/users/ajing/ligandNet/pylib/BindingMoad2011_test/BindingMoad2011/"
__OUTPUT_PROTEINLIGAND__ = "/users/ajing/ligandNet/pylib/proteinChain.txt"

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
        error_obj.write( w.message )
    except:
        error_obj.write( "some value error when reading: " + pdb_code )
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
            if getChainLength( chain ) > 10:
                proteinChains.append( chain.id )
            else:
                for residue in chain.get_list():
                    residueName = residue.id[0]
                    if residueName.startswith("H_"):
                        residueName = residueName[2:]
                    if not residueNameinLigandlist( residueName, ligandlist ):
                        exception_obj.write( "\t".join( [filename, residueName, chain.id] ) + "\n" )
    return proteinChains

def processOneBioUnit( biounitdir ):
    out_obj = open( __OUTPUT_PROTEINLIGAND__, "a")
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL = everyparser.ALL
    structure = makeStructure( biounitdir )
    filename  = biounitdir.split("/")[-1]
    PDBID     = filename[:4].upper()
    PDBChainIDs = returnProteinChainID( structure, filename, __ALL__[PDBID].keys() )
    out_obj.write( filename + "\t" + "".join(PDBChainIDs) + "\n" )
    out_obj.close()

def main():
    # remove error and output files first
    try:
        os.remove( __ERROR_FILE__ )
        os.remove( __EXCEPTION_FILE__ )
        os.remove( __OUTPUT_PROTEINLIGAND__ )
    except:
        pass
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    global __ALL__
    __ALL__ = everyparser.ALL
    pool = Pool(processes = 6)
    argumentlist = [ __BIOUNIT_DIR__ + filename for filename in os.listdir(__BIOUNIT_DIR__)]
    result = pool.map_async( processOneBioUnit, argumentlist)
    resulttxt = result.get()
    print dir(resulttxt)

if __name__ == "__main__":
    main()
