'''
    filter ProbisInput.txt some satisfy aqeel's requirement
'''
from Modular import PROBIS_DIR, PROBIS_COLNAME
from convertToProBiSFormat import assignEachPDBwithNumberofMembers,  assignEachPDBwithNumberofMembers, isleader, ligandFilter
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

Probis_filter = PROBIS_DIR[:-4] + "_filter.txt"

def GetLargeBindingSiteInfo(infile, PDBMemberNumberDict):
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
        ligand  = content[ligand_index]
        pligand = "_".join([ x.split(".")[0] for x in ligand.split("_")])
        if not isleader(PDB, PDBMemberNumberDict):
            continue
        if not PDB in large_dict:
            large_dict[PDB] = dict()
        if not ligand in large_dict[PDB]:
            large_dict[PDB][pligand] = 0
        if bs_size > large_dict[PDB][pligand]:
            large_dict[PDB][pligand] = bs_size
    return large_dict

def LeaderWithLargestBindingSite(infile, PDBMemberNumberDict):
    num_col = len(PROBIS_COLNAME)
    bsdict = GetLargeBindingSiteInfo(infile, PDBMemberNumberDict)
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
        if isleader(PDB, PDBMemberNumberDict):
            if bsdict[PDB][pligand] == bs_size and not onlyOneLigand.checkLigand(PDB, pligand):
                newfile_obj.write(line)
        else:
            newfile_obj.write(line)
    newfile_obj.close()

def main():
    # for everyparser of valid ligand
    everyparser = every_parser("every.csv")
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )
    LeaderWithLargestBindingSite(PROBIS_DIR, pdb_with_numberofmembers)

def test():
    # for everyparser of valid ligand
    everyparser = every_parser("every.csv")
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )
    print GetLargeBindingSiteInfo(PROBIS_DIR, pdb_with_numberofmembers)

if __name__ == "__main__":
    #test()
    main()
