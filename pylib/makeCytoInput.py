"""
    Read probis input and make input file for Cytoscape
"""
from ProbisInputReader import Probis

# addtional info file directory
__ADDITIONAL_INFO__ = "../Data/additionalInfo.txt"
__CYTOSCAPE_INPUT__ = "../Data/CytoscapeInput.txt"
PROTEIN_CLASS = { 1:"OXIDOREDUCTASES", 2:"TRANSFERASES", 3:"HYDROLASES", 4:"LYASES", 5:"ISOMERASES", 6:"LIGASES", 7:"UNCLASSIFIED", 8:"BINDING", 9:"FOLDING", 10:"IMMUNE", 11:"MOBILE", 12:"OTHER", 14:"TOXIN_VIRAL", 15:"TRANSCRIPT_TRANSLATE", 16:"TRANSPORT", 17: "CELL_CYCLE", 18:"SIGNAL_HORMONE", 19:"STRUCTURAL"}
AMINO_ACID = [ "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
NEWNAME_COLUMN = {"INDEX": "BS_ID", "BINDINGSITESIZE": "BS_Size", "ISHEADER": "Leader", "pdbFile": "PDB", "name": "Lig_Name", "affinmicromolar": "Bdata", "ecNumber": "EC_Num"}
NEWNAME_ORDER  = ["BS_ID", "BS_Size", "Leader", "PDB", "EC_Num", "Class", "mw", "Bdata", "Lig_Type", "Lig_Name", "Lig_Size"]

def AdditionalInfoParser( infile ):
    # So here we assume the file has colname
    first_line = 1
    addict = dict()
    for line in open(infile):
        content = [ each.strip() for each in line.strip().split("\t") ]
        if first_line:
            colname = content
            for each in colname:
                addict[each] = []
            first_line = 0
        else:
            for i in range(len(colname)):
                addict[colname[i]].append(content[i])
    return (addict, colname)

def ProteinClass( addict, i ):
    ec_number = addict["ecNumber"][i]
    first_ec = int(ec_number.split(".")[0])
    return PROTEIN_CLASS[first_ec]

def LineWithColumn( dictobj, colname, i ):
    linelist = []
    if colname is None:
        for each in dictobj:
            linelist.append(dictobj[each][i])
    else:
        for each in colname:
            linelist.append(dictobj[each][i])
    linelist = [ str(each) for each in linelist ]
    return "\t".join(linelist)

def LigType(ligandname):
    # take ligand name and return ligand type: Simple, MPart, Peptide
    ligands = ligandname.strip().split()
    if len(ligands) == 1:
        return "Simple"
    if all( i in AMINO_ACID for i in ligands ):
        return "Peptide"
    else:
        return "MPart"

def LigSize(ligandname):
    # return how many different residue in that ligand
    return str(len(ligandname.split()))

def TranslateColname(namelist):
    newlist = []
    for each in namelist:
        if each in NEWNAME_COLUMN:
            newlist.append(NEWNAME_COLUMN[each])
        else:
            newlist.append(each)
    return newlist

def GetOrder(colnames):
    orderlist = []
    for each in colnames:
        orderlist.append(NEWNAME_ORDER.index(each))
    return orderlist

def ReorderList(oldlist, order):
    return [oldlist[i] for i in order]

def MakeCytoInput():
    colname_keep = ["INDEX", "BINDINGSITESIZE", "ISHEADER"]
    probis_dict  = Probis().probisdict
    addit_dir    = __ADDITIONAL_INFO__
    addit_dict, addit_colname   = AdditionalInfoParser(addit_dir)
    addit_colname.pop(0)
    addit_dict["id"] = [int(each) for each in addit_dict["id"]]
    probis_dict["MOADINDEX"] = [int(each) for each in probis_dict["MOADINDEX"]]
    cyto_obj     = open(__CYTOSCAPE_INPUT__, "w")
    ## add first line as column for the file
    all_colnames = colname_keep + addit_colname + ["Class", "Lig_Type", "Lig_Size"]
    all_colnames = TranslateColname(all_colnames)
    colname_order = GetOrder(all_colnames)
    cyto_obj.write( "\t".join(ReorderList(all_colnames, colname_order)) + "\n")
    for i in range(len(probis_dict["MOADINDEX"])):
        moadline  = LineWithColumn( probis_dict, colname_keep, i )
        moadindex = probis_dict["MOADINDEX"][i]
        addit_index = addit_dict["id"].index(moadindex)
        additline = LineWithColumn( addit_dict, addit_colname, addit_index )
        proteinclass = ProteinClass(addit_dict, addit_index)
        ligandname   = addit_dict["name"][addit_index]
        ligand_type  = LigType(ligandname)
        ligand_size  = LigSize(ligandname)
        items        = moadline.split("\t") + additline.split("\t") + [ proteinclass, ligand_type, ligand_size ]
        final_line = "\t".join(ReorderList(items, colname_order)) + "\n"
        cyto_obj.write(final_line)
    cyto_obj.close()


if __name__ == "__main__":
    MakeCytoInput()
