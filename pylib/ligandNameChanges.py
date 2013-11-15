'''
    Process ligand name changes in BindingMOAD
'''

from Modular import LIGAND_NAME_CHANGES
_ligandfile = "LigandNameChanges.txt"

class LigandName(dict):
    def __init__(self):
        self.ligand_container = dict()
        self.ligand_protein_container = dict()

    def __setitem__(self, key, value):
        if len(key) > 1:
            new_key = "_".join(key)
            self.ligand_protein_container[new_key] = value
        else:
            key = key[0]
            self.ligand_container[key] = value

    def __getitem__(self, key, **args):
        if len(key) != 2:
            raise TypeError("must be PDB:ligand to retrieve the value")
        new_key = "_".join(key)
        if new_key in self.ligand_protein_container:
            return self.ligand_protein_container[new_key]
        elif key[1] in self.ligand_container:
            return self.ligand_container[key[1]]
        else:
            return key[1]

def LigandNameFileParser(ligandfile):
    ligand_dict = LigandName()
    for line in open(ligandfile):
        content = line.strip().split()
        value   = content[-1]
        key     = content[:-1]
        ligand_dict[key] = value
    return ligand_dict
_ligand_dict = LigandNameFileParser(_ligandfile)

def GetNewName(PDBID, ligand):
    return _ligand_dict[[PDBID, ligand]]

if __name__ == "__main__":
    print GetNewName("1PHW", "ROB")
    print GetNewName("1PHW", "ICI")
    print GetNewName("1MG5", "NAH")
    print GetNewName("1PHW", "NAH")
