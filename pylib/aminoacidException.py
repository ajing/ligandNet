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

from Modular import EVERYCSV

__AMINO_ACID__ = [ "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

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
    everyparser = every_parser(EVERYCSV)
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
    everyparser = every_parser(EVERYCSV)
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
                # possible ligands
                pos_ligands = [ probisdict["LIGANDID"][x].split(".")[0] for x in indexlist]
                pos_ligands = list(set(pos_ligands))
                if not pos_ligands:
                    pos_ligands = ["not found"]
                if flag:
                    PairCannotFind.append([PDB, ligand] + [","] + pos_ligands)
    print "total number of PDB ligand pairs in every.csv: ", TotalNumberPairs
    print "number of PDB ligand pair cannot find in final result: ", len(PairCannotFind)
    LearningChanges(PairCannotFind)
    pairs_after_filter = PairCannotFindFilter(PairCannotFind)
    print "\n".join([" ".join(each) for each in pairs_after_filter])

def LearningChanges(paircannotfind):
    translate = dict()
    for each in paircannotfind:
        if len(each) == 4:
            ligand = each[1]
            if not each[3] == "not found":
                ch_ligand = each[3]
                translate[ligand] = ch_ligand
    for each in paircannotfind:
        if each[3] == "not found":
            ligand = each[1]
            if ligand in translate:
                newligand = translate[ligand]
                each[3]   = newligand

def NumberofCannotFind(paircannotfindlist):
    count = 0
    for each in paircannotfindlist:
        if each[-1] == "not found":
            count += 1
    print count

def ReLearningCannotFindList(infile):
    pairscannotfind = []
    for line in open(infile):
        before, pos_ligands = line.strip().split(",")
        before = before.strip()
        pos_ligands = pos_ligands.strip()
        beforesplit = before.split()
        PDB = beforesplit[0]
        ligands = beforesplit[1:]
        pairscannotfind.append([PDB, " ".join(ligands)] + [","] + [pos_ligands])
    NumberofCannotFind(pairscannotfind)
    LearningChanges(pairscannotfind)
    NumberofCannotFind(pairscannotfind)
    print "\n".join([" ".join(each) for each in pairscannotfind])

def AddQuestionMark(infile):
    for line in open(infile):
        cont, pos = line.strip().split(",")
        content = cont.split()
        PDB = content[0]
        ligand = content[1:]
        if pos[-1] != "?" and len(ligand) > 1 and not pos[-5:] == "found":
            pos = pos + "?"
        print ",".join([cont, pos])

def ExistingPairParser():
    infile = "tmp"
    pair_list = []
    for line in open(infile):
        PDB, ligand = line.strip().split()
        PDB.strip()
        ligand.strip()
        pair_list.append([PDB, ligand])
    return pair_list

def PairCannotFindFilter(paircannotfind):
    newlist = []
    existpairs = ExistingPairParser()
    print existpairs
    print len(existpairs)
    for each in paircannotfind:
        ligand = each[1]
        if each in existpairs:
            continue
        #if len(ligand.split()) > 1:
        #    continue
        #if len(ligand) > 2 and ligand[2] == "P" and ligand[0] == "A":
        #    continue
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
    print "\n".join([" ".join(each) for each in pairs_after_filter])

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
    ProbisInputStat()
    # Help to increase the speed of filtering pair cannot find file
    #ReLearningCannotFindList("outfile4")
    #SimpleFilter("outfile6")
    #AddQuestionMark("outfile6")
