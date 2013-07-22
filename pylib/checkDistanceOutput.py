"""
    This program is to validate contact residue for ligand from Rich's code amino_acid_counter.pl
    Here chain is ligand, atom is the atom on protein pocket
    Author: ajing
    Date:   7/20/2013
"""
from compareMOADBiounit import makeStructure
import numpy

# exception column which is wrong
__WRONG_COLUMN__ = "/users/ajing/ligandNet/pylib/wrongRichsColumn.txt"

# The original colname of file
colname = [
    "BIOUNIT",
    "BIOUNITFILE",
    "ligandName",
    "ligandChainID",
    "ligandChainIDNEW",
    "proteinChainID",
    "residueName",
    "residueNumber",
    "insertion",
    "AtomName",
    "AtomNumber",
    "Distance",
    "Metal"
]


def processEachLine( filename, colname ):
    wrong_obj = open( __WRONG_COLUMN__, 'w' )
    for eachline in open( filename ):
        row_dict = dict()
        content  = eachline.split(",")
        content_len = len( content )
        col_len     = len( colname )
        if content_len != col_len:
            raise "corrupted row:" + eachline
        for i in range( col_len ):
            row_dict[ colname[i] ] = content[i]
        biounit        = row_dict[ "BIOUNITFILE" ]
        structure      = makeStructure( biounit )
        ligandChainID  = row_dict[ "ligandChainID" ]
       # print biounit
       # print eachline
       # print ligandChainID
        if not ligandChainID:
            continue
        chain          = getChain( structure, ligandChainID )
        atomChainID    = row_dict[ "proteinChainID" ]
        atomResidue    = int( row_dict[ "residueNumber" ])
        atomchain      = getChain( structure, atomChainID )
        atomName       = row_dict[ "AtomName" ]
        try:
            atom       = getAtom( structure, atomChainID, atomResidue, atomName )
        except:
            atomResidueName = row_dict[ "residueName" ]
            atomResidueFull = ( "H_" + atomResidueName.rjust(3), atomResidue, ' ')
           # print atomchain.child_dict
           # print atomResidueFull
            atom       = getAtom( structure, atomChainID, atomResidueFull, atomName )
        if minimunDistanceChainToAtom( chain, atom ) > 4:
            wrong_obj.write( eachline )
    wrong_obj.close()

def getChain( structure, chainID):
    # here chain is ligand
    # choose model 0 with chainID
    for eachmodel in structure:
        try:
            return eachmodel[chainID]
        except:
            continue
    raise "cannot find chain"

def getAtom( structure, chainID, residueNumber, atomName ):
    # return the atom object in that structure
    for eachmodel in structure:
        try:
            return eachmodel[chainID][residueNumber][atomName]
        except:
            continue
    raise "cannot find atom"

def displayChildren( anode ):
    for each in anode.get_list():
        print each.get_full_id()
        print each.id

def dist(x,y):
    return numpy.sqrt(numpy.sum((x-y)**2))

def minimunDistanceChainToAtom( Chain, Atom ):
    minimun = 100
    for eachResidue in Chain.get_list():
        for eachAtom in eachResidue.get_list():
            distance = dist( eachAtom.coord, Atom.coord )
            if minimun > distance:
                minimun = distance
    return minimun


def debug():
    biounit = "/users/ajing/ligandNet/pylib/BindingMoad2011_test/BindingMoad2011/3l8h.bio1"
    structure = makeStructure( biounit )
    chainid   = "H"
    chain     = getChain( structure, chainid )
    atomchainid = "A"
    atomresidue = 53
    atomname = "OG"
    atom = getAtom( structure, atomchainid, atomresidue, atomname )
    print atom
    print atom.coord
    print minimunDistanceChainToAtom( chain, atom )

if __name__ == "__main__":
    richfile = "/users/ajing/ligandNet/pylib/tmp_test/final.txt"
    processEachLine( richfile, colname )
