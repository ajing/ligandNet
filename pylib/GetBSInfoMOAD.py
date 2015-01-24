'''
    This file is for convert original file of Rich's amino_acid_counter.pl to only valid ligand binding sites
    author: ajing
    Date: 1/23/2015
'''

import sys
import csv
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from Modular import PROBIS_DIR, BIOUNIT_DIR, RICHOUT_DIR, PROTEIN_DIR, MOADINDEX, OBSOLETE, EVERYCSV, RICH_COLNAME
from every_parser import every_parser
from convertToProBiSFormat import ligandFilter, processStrangeLigandName, checkValid
from RichOutParser import RichOutParser

__OUTPUT__ = "bs_info.txt"

alpha = 0.5   # only when ligands contains more than certain percent as the ligand name, I consider this ligand can map to certain ligand in MOAD


class GetOtherRichData:
    def __init__(self):
        allrich = dict()
        for eachitem in csv.reader(open(RICHOUT_DIR)):
            pdbid = eachitem[RICH_COLNAME.index("BIOUNIT")]
            if pdbid in allrich:
                allrich[pdbid].append(eachitem)
            else:
                allrich[pdbid] = [eachitem]
        self.allrich = allrich

    def getResName(self, biounit, pro_chainid, resnum):
        for each in self.allrich[biounit]:
            if each[RICH_COLNAME.index("proteinChainID")] == pro_chainid and int(each[RICH_COLNAME.index("residueNumber")]) == resnum:
                return each[RICH_COLNAME.index("residueName")]

    def getDist(self, biounit, pro_chainid, resnum):
        mindist = 1000
        for each in self.allrich[biounit]:
            if each[RICH_COLNAME.index("proteinChainID")] == pro_chainid and int(each[RICH_COLNAME.index("residueNumber")]) == resnum:
                dist = float(each[RICH_COLNAME.index("Distance")])
                if dist < mindist:
                    mindist = dist
        return mindist


# Output format: BioUnitID, ligandChainID, proteinChainID, residueNumber
def makeBSInput( ProBiS_dict, validliganddict, outfile):
    # Basically in the following format:
    # Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
    # Output2 format: Index, ligandName.ligandChainID,.. ( no info for ligand residue, so ignore that )
    # Example for output2: 00001 ATP.123.A,ALA.43.C,TYR.44.C
    # 7/30/2013 only keep one ligand for one PDB
    out_obj = open( outfile, "w" )
    # 7/30/2013 Here is the structure to keep track only one case for each ligand, remove the duplication in biounit also
    onlyOneLigandEntry = ligandFilter()
    richdata = GetOtherRichData()
    for eachBioUnit in sorted( ProBiS_dict.keys() ):
        BioUnitID   = eachBioUnit
        PDBID       = BioUnitID.split('.')[0]
        # 11/12/2013 for ignoring obsolete protein
        if PDBID in OBSOLETE:
            continue
        for eachligand in sorted(ProBiS_dict[eachBioUnit].keys()):
            # 9/6 sort protein chain id also
            ligandName, ligandChain = eachligand.split('.')
            for eachproteinChainID in sorted( ProBiS_dict[eachBioUnit][eachligand].keys() ):
                ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID].sort()
                for resnum in ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID]:
            # 7/30/2013 for only one ligand for leader
            #if isleader(PDBID, PDBMemberNumberDict) and onlyOneLigandEntry.checkLigand( PDBID, ligandName ):
            #    continue
                    if checkValid( PDBID, ligandName, validliganddict ):
                        resnam = richdata.getResName(BioUnitID, eachproteinChainID, resnum)
                        resdist = richdata.getDist(BioUnitID, eachproteinChainID, resnum)
                        if resnam is None or resdist == 1000:
                            print BioUnitID, ligandName, eachproteinChainID, str(resnum)
                        content = [PDBID, BioUnitID, ligandName, eachproteinChainID, resnam, resnum, resdist]
                        out_obj.write("\t".join(map(str, content)) + "\n")
    out_obj.close()

def main():
    everyparser = every_parser(EVERYCSV)
    everyparser.find_PDBID_ValidLigand()
    infiledir = RICHOUT_DIR
    richout     = RichOutParser( infiledir, everyparser.ALL )
    probisdict  = richout.obj
    makeBSInput( probisdict, everyparser.ALL, __OUTPUT__)

def test():
    other = GetOtherRichData()
    print other.getResName("4K3N.BIO2", "L1", 3)
    print other.getDist("4K3N.BIO2", "L1", 3)
    print other.getDist("1A9X.BIO1", "E5", 911)
    #print other.getResName("10GS.BIO1", "A", 10)
    #print other.getResName("10GS.BIO1", "A", 13)
    #print other.getResName("10GS.BIO1", "A", 38)
    #print other.getResName("9XIM.BIO1", "D", 396)


if __name__ == "__main__":
    main()
    #test()
