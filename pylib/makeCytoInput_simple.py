"""
    Read probis input and make input file for Cytoscape
"""
from ProbisInputReader import Probis
from makeCytoInput import LineWithColumn, LigType, LigSize, PROTEIN_CLASS
__CYTOSCAPE_INPUT__ = "../Data/CytoscapeInput.txt"

## For every_parser
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

def RefineLigand(ligandname):
    return " ".join([ x.split(".")[0] for x in ligandname.split("_") ])

def MakeCytoInput():
    colname = ["BS_ID", "BS_Size", "PDB", "EC_Num", "Class", "Lig_Type", "Lig_Name", "Lig_size"]
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    probis_dict  = Probis().probisdict
    cyto_obj     = open(__CYTOSCAPE_INPUT__, "w")
    cyto_obj.write("\t".join(colname) + "\n")
    for i in range(len(probis_dict["MOADINDEX"])):
        index   = probis_dict["INDEX"][i]
        biounit = probis_dict["BIOUNIT"][i]
        pdb     = biounit.split(".")[0].upper()
        ligand  = probis_dict["LIGANDID"][i]
        ligand  = RefineLigand(ligand)
        ligtype = LigType(ligand)
        bs_size = probis_dict["BINDINGSITESIZE"][i]
        num_member = probis_dict["HEADER"][i]
        lig_size= LigSize(ligand)
        ec_number= everyparser.ALL_class[pdb]
        first_ec = int(ec_number.split(".")[0])
        pro_class= PROTEIN_CLASS[first_ec]
        line = [index, bs_size, pdb, ec_number, pro_class, ligtype, ligand, lig_size]
        cyto_obj.write("\t".join(line) + "\n")
    cyto_obj.close()

if __name__ == "__main__":
    MakeCytoInput()
