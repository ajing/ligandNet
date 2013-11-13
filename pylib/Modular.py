"""
    Data  Warehouse
"""
#PROBIS_DIR = "./ProbisInput_aqeel.txt"
PROBIS_DIR = "./ProbisInput_test.txt"
BIOUNIT_DIR = "../2012_biounits/"
RICHOUT_DIR = "../tmp2012/final.txt"
PROTEIN_DIR = "../Data/proteinChain_2012.txt"
PROBIS_COLNAME = [
                "INDEX",
                "MOADINDEX",
                "BIOUNIT",
                "LIGANDID",
                "PROTEINCHAIN",
                "BINDINGSITE",
                "BINDINGSITESIZE",
                "HEADER",
                "ISHEADER"
                ]

# BindingMOAD numbering for PDB and binded ligands
MOADINDEX = "../Data/MOADindex.txt"

# OBSOLETE LIST
OBSOLETE = ['1UVW', '2BPL', '2EBC', '2G8O', '2OFI', '2ZCD','2ZCP', '3C67', '3C68', '3C69', '3N3V', '3QD1', '3RBO', '3TVS', '3U7F', '3U7G', '3U7H', '3UG9', '3ZVM', '3ZVN', '4A75', '4A76' ]
