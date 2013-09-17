"""
    This code is to deal with special amino acid ligands which is hard to distinguish the distance.
    Author: ajing
    Date  : 9/9/2013

"""

from ProbisInputReader import Probis
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

__AMINO_ACID__ = [ "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

def numberOfAminoAcidLigand():
    probis   = Probis()
    ligandID = probis.probisdict["LIGANDID"]
    n = 0
    for index in range(len(ligandID)):
        if "_" in ligandID[index]:
            continue
        ligandName = ligandID[index].split(".")[0]
        if ligandName in __AMINO_ACID__:
            n = n + 1
            print probis.probisdict["BIOUNIT"][index] + "," + ligandID[index] + "," + probis.probisdict["BINDINGSITESIZE"][index]
    print n
    #indexlist = []
    #for each in probis.probisdict["MOADINDEX"]:
    #    if not each in indexlist:
    #        indexlist.append( int(each) )
    #print max( indexlist )
    #print len( indexlist )
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL         = everyparser.ALL
    validligand = []
    for PDB in ALL:
        for ligand in ALL[PDB]:
            if ALL[PDB][ligand]:
                validligand.append( ligand )
    print validligand
    print len(validligand)

def ligandNotInFile():
    probis   = Probis()
    PDBIDs   = [ each.split(".")[0].upper() for each in probis.probisdict["BIOUNIT"] ]
    PDBIDs_len = len(PDBIDs)
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL         = everyparser.ALL
    longligandCount = 0
    PairCannotFind  = []
    for PDB in ALL:
        indexlist = [ each for each in range(PDBIDs_len) if PDB == PDBIDs[each] ]
        for ligand in ALL[PDB]:
            # flag for whether such pair exists in file
            flag = 1
            if ALL[PDB][ligand]:
                for index in indexlist:
                    ligandID = probis.probisdict["LIGANDID"][index]
                    ligandName = ligandID.split(".")[0]
                    if ligandName in ligand or ligand in ligandName:
                        flag = 0
                        break
                if flag:
                    PairCannotFind.append([PDB, ligand])
    print longligandCount
    print PairCannotFind
    print len(PairCannotFind)

if __name__ == "__main__":
    #numberOfAminoAcidLigand()
    ligandNotInFile()

