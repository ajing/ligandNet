"""
    Convert to Probis input, only keep all leader ligands
"""

from convertToProBiSFormat import *
from ProbisInputReader import Probis
from Modular import PROBIS_DIR

def getExistID(infile, astring):
    for line in open(infile):
        index = line.strip().split("\t")[0]
        if astring in line:
            return int(index)

def makeProBiSInputLeaderOnly( ProBiS_dict, validliganddict, outfile, previousfile, outfile2, BioChainsDict = None, NumberDict = None, PDBMemberNumberDict = None ):
    # Basically in the following format:
    # Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
    # Output2 format: Index, ligandName.ligandChainID,.. ( no info for ligand residue, so ignore that )
    # Example for output2: 00001 ATP.123.A,ALA.43.C,TYR.44.C
    # 10/25
    out_obj = open( outfile, "a" )
    ExistingPairs = []
    # 7/30/2013 Here is the structure to keep track only one case for each ligand, remove the duplication in biounit also
    onlyOneLigandEntry = ligandFilter()
    # 9/17 read existing output and generate continuous index
    existingContent, index = readExistingOutput(previousfile)
    indexG = indexGenerator( 75000 )
    for eachBioUnit in sorted( ProBiS_dict.keys() ):
        BioUnitID   = eachBioUnit
        PDBID       = BioUnitID.split('.')[0]
        for eachligand in sorted(ProBiS_dict[eachBioUnit].keys()):
            oneLine       = oneLineInfo()
            # The total number of binding site residue
            bindingSiteNumber = 0
            # 9/6 sort protein chain id also
            for eachproteinChainID in sorted( ProBiS_dict[eachBioUnit][eachligand].keys() ):
                ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID].sort()
                bindingSiteNumber = bindingSiteNumber + len(ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID])
                chainInfo = eachproteinChainID + " and (" + ",".join( map( str, ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID] ) ) + ")"
                oneLine.OneMoreChain( chainInfo )
            oneLine.OneMoreChain()
            ligandName, ligandChain = eachligand.split('.')
            # 7/30/2013 for only one ligand
            if not isheader(PDBID, PDBMemberNumberDict):
                continue
            if not BioChainsDict.has_key(BioUnitID.lower()):
                continue
            if checkValid( PDBID, ligandName, validliganddict ):
                ligandChainID = processStrangeLigandName( eachligand )
                if BioChainsDict is None:
                    # whether include binding PDB chains or not
                    oneLine.addLeft("\t".join( [ BioUnitID.lower(), ligandChainID, allproteinChainInfo ] ) )
                else:
                    numbering = getBindingMoadID( PDBID.upper(), ligandName, NumberDict )
                    oneLine.addLeft("\t".join( [ numbering, BioUnitID.lower(), ligandChainID, BioChainsDict[BioUnitID.lower()] ]) )
                    oneLine.addRight("\t".join( [ str(bindingSiteNumber), PDBMemberNumberDict[PDBID.upper()] ])  )
                if not oneLine.string in existingContent:
                    oneLine.setIndex( indexG.next() )
                    oneLine.writeToFile( out_obj )
                else:
                    oneLine.setIndex( getExistID(previousfile, oneLine.string) )
                    oneLine.writeToFile( out_obj )
            # add one index for each binding site/ligand
    out_obj.close()

def main():
    # for everyparser of valid ligand
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )

    infiledir = "/users/ajing/ligandNet/tmp_test/final.txt"
    outfiledir = "/users/ajing/ligandNet/pylib/ProbisInputLeader.txt"
    outfiledir2 = "/users/ajing/ligandNet/Data/IndexwithLigandInfo.txt"
    ## file for protein id : protein chains
    proteinChainFile = "/users/ajing/ligandNet/Data/proteinChain.txt"
    proteinchaindict = returnChainsForPDBID( proteinChainFile )
    ## file for aqeel's numbering for protein ligand pair
    aqeelfile = "/users/ahmedaq/work/Probis/LigandID_PDB_LigName_tab.nosql"
    aqeeldict   = aqeelNumberingParse( aqeelfile )
    richout     = RichOutParser( infiledir, everyparser.ALL )
    probisdict  = richout.obj
    makeProBiSInputLeaderOnly( probisdict, everyparser.ALL, outfiledir, PROBIS_DIR, outfiledir2, proteinchaindict, aqeeldict, pdb_with_numberofmembers)

if __name__ == "__main__":
    main()
