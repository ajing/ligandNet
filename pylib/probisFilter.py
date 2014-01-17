'''
    filter ProbisInput.txt some satisfy aqeel's requirement
'''
from Modular import PROBIS_DIR, PROBIS_COLNAME
from convertToProBiSFormat import assignEachPDBwithNumberofMembers,  assignEachPDBwithNumberofMembers, isleader, ligandFilter, checkValid
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

Probis_filter = PROBIS_DIR[:-4] + "_filter.txt"

def GetLargeBindingSiteInfo(infile, PDBMemberNumberDict, validdict):
    num_col = len(PROBIS_COLNAME)
    bs_size_index = PROBIS_COLNAME.index("BINDINGSITESIZE")
    biounit_index = PROBIS_COLNAME.index("BIOUNIT")
    ligand_index  = PROBIS_COLNAME.index("LIGANDID")
    large_dict    = dict()
    for line in open(infile):
        content = line.strip().split("\t")
        if len(content) != num_col:
            raise TypeError("dismatch of column number:" + line)
        bs_size = int(content[bs_size_index])
        biounit = content[biounit_index]
        PDB     = biounit.split(".")[0].upper()
        ligand  = content[ligand_index]  # ligand from ProbisInput.txt LIG.A_LAG.B
        pligand = "_".join([ x.split(".")[0] for x in ligand.split("_")]) # pure ligand name LIG_LAG
        if not isleader(PDB, PDBMemberNumberDict):
            continue
        if not PDB in large_dict:
            large_dict[PDB] = dict()
        formal_ligand = checkValid(PDB, " ".join(pligand.split("_")), validdict)
        if not formal_ligand in large_dict[PDB]:
            large_dict[PDB][formal_ligand] = 0
        if bs_size > large_dict[PDB][formal_ligand]:
            large_dict[PDB][formal_ligand] = bs_size
    return large_dict

def LeaderWithLargestBindingSite(infile, PDBMemberNumberDict, validdict):
    num_col = len(PROBIS_COLNAME)
    bsdict = GetLargeBindingSiteInfo(infile, PDBMemberNumberDict, validdict)
    bs_size_index = PROBIS_COLNAME.index("BINDINGSITESIZE")
    biounit_index = PROBIS_COLNAME.index("BIOUNIT")
    ligand_index  = PROBIS_COLNAME.index("LIGANDID")
    onlyOneLigand = ligandFilter()
    newfile_obj   = open(Probis_filter, "w")
    for line in open(infile):
        content = line.strip().split("\t")
        if len(content) != num_col:
            raise TypeError("dismatch of column number:" + line)
        bs_size = int(content[bs_size_index])
        biounit = content[biounit_index]
        PDB     = biounit.split(".")[0].upper()
        ligand  = content[ligand_index]
        pligand = "_".join([ x.split(".")[0] for x in ligand.split("_")])
        formal_ligand = checkValid(PDB, " ".join(pligand.split("_")), validdict)
        if isleader(PDB, PDBMemberNumberDict):
            if bsdict[PDB][formal_ligand] == bs_size and not onlyOneLigand.checkLigand(PDB, pligand):
                newfile_obj.write(line)
        else:
            newfile_obj.write(line)
    newfile_obj.close()

def main():
    # for everyparser of valid ligand
    everyparser = every_parser("every.csv")
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )
    LeaderWithLargestBindingSite(PROBIS_DIR, pdb_with_numberofmembers, everyparser.ALL)

def test():
    # for everyparser of valid ligand
    everyparser = every_parser("every.csv")
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )
    print GetLargeBindingSiteInfo(PROBIS_DIR, pdb_with_numberofmembers)

if __name__ == "__main__":
    #test()
    main()
