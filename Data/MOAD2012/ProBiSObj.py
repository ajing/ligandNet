'''
    Create an object to contain all data for each binding site
'''

from convertToProBiSFormat import oneLineInfo

class BindingSite:
    def __init__(self, biounit, ligand, ligandchains):
        self.id      = 0
        self.biounit = biounit
        self.ligand  = ligand
        self.ligandchains  = ligandchains
        self.proteinchains = proteinchains
        self.PDB = biounit.split(".")[0].upper()
        self.BS_size = 0
    pass

def CreateBSObj(Probis_dict):
    BS_list = []
    for eachBioUnit in sorted(Probis_dict.keys()):
        for eachligand in sorted(Probis_dict[eachBioUnit].keys()):
            bindingsite_size = 0
            ligandName, ligandChain = eachligand.split('.')
            one_BS = BindingSite(eachBioUnit, ligandName, ligandChain)
            for eachproteinChainID in sorted( ProBiS_dict[eachBioUnit][eachligand].keys() ):

