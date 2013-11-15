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

# the pdb we cannot find any processed
__pdbcannotfind__  = ""

from RichOutParser import RichOutParser
from Modular import BIOUNIT_DIR

def numberOfAminoAcidLigand():
    probis   = Probis()
    ligandID = probis.probisdict["LIGANDID"]
    aminoacid = dict()
    for index in range(len(ligandID)):
        if "_" in ligandID[index]:
            continue
        ligandName = ligandID[index].split(".")[0]
        if ligandName in __AMINO_ACID__:
            PDB = probis.probisdict["BIOUNIT"][index].split(".")[0].upper()
            ligand = ligandID[index].split(".")[0]
            try:
                if not ligand in aminoacid[PDB]:
                    aminoacid[PDB].append(ligand)
            except:
                aminoacid[PDB] = [ligand]
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

def ligandNotInFile(probisdict):
    PDBIDs   = [ each.split(".")[0].upper() for each in probisdict["BIOUNIT"] ]
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
                    ligandID = probisdict["LIGANDID"][index]
                    ligandName = ligandID.split(".")[0]
                    if ligandName in ligand or ligand in ligandName:
                        flag = 0
                        break
                if flag:
                    PairCannotFind.append([PDB, ligand])
    print "total number of PDB ligand pairs in every.csv: ", TotalNumberPairs
    print "number of PDB ligand pair cannot find in final result: ", len(PairCannotFind)

def PairCannotFindFilter(paircannotfind):
    newlist = []
    for each in paircannotfind:
        ligand = each[1]
        if len(ligand.split()) > 1:
            continue
        if len(ligand) > 2 and ligand[2] == "P" and ligand[0] == "A":
            continue
        newlist.append(each)
    return newlist

def ExtendLigandNameChanges(pairscannotfind):
    infile = "LigandNameChanges.txt"
    ligandcannot_list = []
    for each in open(infile):
        content = each.strip().split()
        if len(content) > 2:
            ligandchanged = content[1]
            for eachcannot in pairscannotfind:
                if eachcannot[1] == ligandchanged:
                    ligandcannot_list.append(eachcannot)
    return sorted(ligandcannot_list, key = lambda k: k[1])

def findBioUnitfromPDB(biounitlist, PDB):
    biounit_match = []
    for eachbiounit in biounitlist:
        if eachbiounit[:4].upper() == PDB:
            biounit_match.append(eachbiounit)
    return biounit_match

def ligandExistInList(aligand, ligandlist):
    for each in ligandlist:
        if each in aligand or aligand in each:
            return True
    return False

def ligandNotInFileRichOut(richdict):
    biounits = richdict.keys()
    PDBIDs   = [ each.split(".")[0].upper() for each in richdict.keys() ]
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ALL         = everyparser.ALL
    PairCannotFind  = []
    for PDB in ALL:
        for ligand in ALL[PDB]:
            flag = 1
            # flag for whether such pair exists in file
            if ALL[PDB][ligand]:
                biounit_list  = findBioUnitfromPDB(biounits, PDB)
                for eachbiounit in biounit_list:
                    ligand_names = [x.split(".")[0] for x in richdict[eachbiounit].keys()]
                    if ligandExistInList(ligand, ligand_names):
                        flag = 0
                        break
                if flag:
                    PairCannotFind.append([PDB, ligand])
    print "number of PDB ligand pair cannot find in final result: ", len(PairCannotFind)
    pairs_after_filter = PairCannotFindFilter(PairCannotFind)
    print ExtendLigandNameChanges(pairs_after_filter)

def RichOutStatistics( richoutdict ):
    import os
    biounits = richoutdict.keys()
    print "Number of biounit files: " + str(len(biounits))
    PDBs     = set( [ each.split(".")[0] for each in biounits ])
    print "Number of PDBs: " + str(len(PDBs))
    filelist = os.listdir( BIOUNIT_DIR )
    print "Total number of biounit files: " + str(len(filelist))
    AllPDBs     = set( [ each.split(".")[0] for each in filelist])
    print "Total number of PDBs: " + str(len(AllPDBs))
    ligandNotInFileRichOut(richoutdict)

def ProbisInputStat():
    import os
    listfiles = os.listdir(BIOUNIT_DIR)
    probis   = Probis()
    biounits = probis.probisdict["BIOUNIT"]
    allbiounitfiles = listfiles
    biounits = list(set(biounits))
    print "Number of biounit files in ProbisInput.txt: " + str(len(biounits))
    PDBs     = set( [ each.split(".")[0] for each in biounits ])
    print "Number of PDBs: " + str(len(PDBs))
    ligandNotInFile(probis.probisdict)


if __name__ == "__main__":
    #numberOfAminoAcidLigand()
    richdict = RichOutParser("final.txt")
    richdict = richdict.obj
    RichOutStatistics( richdict )
    # probis input stat
    #ProbisInputStat()
