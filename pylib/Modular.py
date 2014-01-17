"""
    Data  Warehouse
    Please call all files in corresponding data directory
"""
PROBIS_DIR = "ProbisInput.txt"
#PROBIS_DIR = "./ProbisInput_test.txt"
BIOUNIT_DIR = "../../2012_biounits/"
#BIOUNIT_DIR = "../../BindingMoad2011_test/BindingMoad2011/"
RICHOUT_DIR = "final.txt"
PROTEIN_DIR = "proteinChain.txt"
PROBIS_COLNAME = [
                "INDEX",
                "MOADINDEX",
                "BIOUNIT",
                "LIGANDID",
                "PROTEINCHAIN",
                "BINDINGSITE",
                "BINDINGSITESIZE",
                "HEADER",
                "ISHEADER",
# This is for MOAD2012
                "BINDINGCHAIN"
                ]

# BindingMOAD numbering for PDB and binded ligands
MOADINDEX = "MOADindex.txt"

# OBSOLETE LIST
OBSOLETE = ['1UVW', '2BPL', '2EBC', '2G8O', '2OFI', '2ZCD','2ZCP', '3C67', '3C68', '3C69', '3N3V', '3QD1', '3RBO', '3TVS', '3U7F', '3U7G', '3U7H', '3UG9', '3ZVM', '3ZVN', '4A75', '4A76' ]

# Ligand name changes file
# this file is obsolte after Rich correct ligand names for MOAD2012
#LIGAND_NAME_CHANGES = "LigandNameChanges.txt"

# every.csv
EVERYCSV = "every.csv"
