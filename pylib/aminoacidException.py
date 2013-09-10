"""
    This code is to deal with special amino acid ligands which is hard to distinguish the distance.
    Author: ajing
    Date  : 9/9/2013

"""

from ProbisInputReader import Probis

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

if __name__ == "__main__":
    numberOfAminoAcidLigand()
