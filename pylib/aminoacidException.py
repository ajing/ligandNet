"""
    This code is to deal with special amino acid ligands which is hard to distinguish the distance.
    Basically the ligands which cannot be found by Richard's code
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
    aminoacid = dict()
    for index in range(len(ligandID)):
        if "_" in ligandID[index]:
            continue
        ligandName = ligandID[index].split(".")[0]
        if ligandName in __AMINO_ACID__:
            #print probis.probisdict["BIOUNIT"][index] + "," + ligandID[index] + "," + probis.probisdict["BINDINGSITESIZE"][index]
            PDB = probis.probisdict["BIOUNIT"][index].split(".")[0].upper()
            ligand = ligandID[index].split(".")[0]
            try:
                if not ligand in aminoacid[PDB]:
                    aminoacid[PDB].append(ligand)
            except:
                aminoacid[PDB] = [ligand]
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
    aminoligandinProBiS = []
    for PDB in ALL:
        for ligand in ALL[PDB]:
            if ALL[PDB][ligand] and ligand in __AMINO_ACID__:
                validligand.append( ligand )
            if PDB in aminoacid and ligand in aminoacid[PDB]:
                aminoligandinProBiS.append(PDB + ligand)
    print "total number of peptide ligand(amino acide as component) in ProbisInput: ", len(aminoligandinProBiS)
    print "number of valid peptide:", len(validligand)

def ligandNotInFile():
    probis   = Probis()
    PDBIDs   = [ each.split(".")[0].upper() for each in probis.probisdict["BIOUNIT"] ]
    PDBIDs_len = len(PDBIDs)
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL         = everyparser.ALL
    PairCannotFind  = []
    TotalNumberPairs = 0
    for PDB in ALL:
        indexlist = [ each for each in range(PDBIDs_len) if PDB == PDBIDs[each] ]
        for ligand in ALL[PDB]:
            # flag for whether such pair exists in file
            flag = 1
            if ALL[PDB][ligand]:
                TotalNumberPairs = TotalNumberPairs + 1
                for index in indexlist:
                    ligandID = probis.probisdict["LIGANDID"][index]
                    ligandName = ligandID.split(".")[0]
                    if ligandName in ligand or ligand in ligandName:
                        flag = 0
                        break
                if flag:
                    PairCannotFind.append([PDB, ligand])
    print "total number of PDB ligand pairs in every.csv: ", TotalNumberPairs
    print "number of PDB ligand pair cannot find in final result: ", len(PairCannotFind)
    print PairCannotFind

if __name__ == "__main__":
    #numberOfAminoAcidLigand()
    ligandNotInFile()
