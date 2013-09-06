'''
    This file is for convert original file of Rich's amino_acid_counter.pl to what ProBiS like
    author: ajing
    Data  : 7/15/2013
    ChangeLog:
        7/29/2013 add another column for number of members in that family
'''

## For every_parser
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

# The original colname of file
colname = [
    "BIOUNIT",
    "BIOUNITFILE",
    "ligandName",
    "ligandChainID",
    "ligandChainIDNEW",
    "proteinChainID",
    "residueName",
    "residueNumber",
    "insertion",
    "AtomName",
    "AtomNumber",
    "Distance",
    "Metal"
]

# Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
# like "[:A and (12, 13, 14, 15)]"

def convertProBiS( infile ):
    output_dict = dict()
    for line in open(infile):
        content = line.strip().split(",")
        content_dict = dict()
        tot_num = len( content )
        if tot_num != len( colname ):
            raise "Some problem with this line: " + content
        for idx in range( tot_num ):
            content_dict[ colname[idx] ] = content[ idx ]
        if content_dict[ "proteinChainID" ] == "z":
            continue
        try:
            output_dict[ content_dict["BIOUNIT"] ]
        except:
            output_dict[ content_dict["BIOUNIT"] ] = dict()
        try:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ]
        except:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ] = dict()
        try:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ]
        except:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ]  = []
        if not int(content_dict["residueNumber"]) in output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ]:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ].append( int(content_dict["residueNumber"]) )
    return output_dict

def processStrangeLigandName( completeLigandName ):
    # Here completeLigandName means like BGC.D or GLC BGC.B_C
    ligandName = completeLigandName.split(".")[0]
    ChainID    = completeLigandName.split(".")[1]
    ligandName.strip()
    ligandNamelist = ligandName.split(" ")
    ChainIDlist    = ChainID.split("_")
    ligandNamelen  = len(ligandNamelist)
    if ligandNamelen == 1:
        return completeLigandName
    else:
        newNamelist = []
        for i in range(ligandNamelen):
            try:
                newNamelist.append( ligandNamelist[i] + "." + ChainIDlist[i] )
            except:
                raise "Strange ligand name: " + completeLigandName
        return "_".join( newNamelist )

def returnChainsForPDBID( BioUnitChainsDIR ):
    # Should return a dictionary object with biounit file name as key and all chains in that biounit file as value
    BioUnitChainDict = dict()
    for line in open( BioUnitChainsDIR ):
        content = line.strip().split( "\t" )
        try:
            Biounitfile = content[0].lower()
            Chains      = content[1]
            BioUnitChainDict[Biounitfile] = Chains
        except:
            print "cannot find protein chains for " + content[0]
    return BioUnitChainDict

def makeProBiSInput( ProBiS_dict, validliganddict, outfile, outfile2, BioChainsDict = None, NumberDict = None, PDBMemberNumberDict = None ):
    # Basically in the following format:
    # Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
    # Output2 format: Index, ligandName.ligandChainID,.. ( no info for ligand residue, so ignore that )
    # Example for output2: 00001 ATP.123.A,ALA.43.C,TYR.44.C
    # 7/30/2013 only keep one ligand for one PDB
    out_obj = open( outfile, "w" )
    index = 1
    ExistingPairs = []
    # 7/30/2013 Here is the structure to keep track only one case for each ligand, remove the duplication in biounit also
    ligandtrack = dict()
    for eachBioUnit in ProBiS_dict.keys():
        ligandlist  = ",".join(ProBiS_dict[eachBioUnit].keys())
        BioUnitID   = eachBioUnit
        PDBID       = BioUnitID.split('.')[0]
        if not PDBID in ligandtrack.keys():
            ligandtrack[PDBID] = []
        try:
            ligandList = ProBiS_dict[eachBioUnit].keys()
            ligandList.sort()
            print ligandList
        except:
            print "something wrong"
            continue
        for eachligand in ligandList:
            Index = str(index).zfill(5)
            ligandChainID = eachligand
            allproteinChainInfo = ""
            # The total number of binding site residue
            bindingSiteNumber = 0
            for eachproteinChainID in ProBiS_dict[eachBioUnit][eachligand].keys():
                ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID].sort()
                bindingSiteNumber = bindingSiteNumber + len(ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID])
                if allproteinChainInfo:
                    allproteinChainInfo = allproteinChainInfo + " or :" + eachproteinChainID + " and (" + ",".join( map( str, ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID] ) ) + ")"
                else:
                    allproteinChainInfo = allproteinChainInfo + "[:" + eachproteinChainID + " and (" + ",".join( map( str, ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID] ) ) + ")"
            allproteinChainInfo = allproteinChainInfo + "]"
            ligandName  = ligandChainID.split('.')[0]
            ligandChain = ligandChainID.split('.')[1]
            # 7/30/2013 for only one ligand
            if ligandName in ligandtrack[PDBID]:
                continue
            else:
                ligandtrack[PDBID].append(ligandName)
            try:
                if validliganddict[PDBID.upper()][ligandName]:
                    ligandChainID = processStrangeLigandName( ligandChainID )
                    if BioChainsDict is None:
                        # whether include binding PDB chains or not
                        aline = "\t".join( [ Index, BioUnitID.lower(), ligandChainID, allproteinChainInfo ] ) + "\n"
                        index = index + 1
                    else:
                        try:
                            NumberingKey = "\t".join([PDBID.upper(), ligandName])
                            numbering = NumberDict[NumberingKey]
                            aline = "\t".join( [ Index, numbering, BioUnitID.lower(), ligandChainID, BioChainsDict[ BioUnitID.lower() ], allproteinChainInfo, str(bindingSiteNumber), PDBMemberNumberDict[PDBID.upper()] ] ) + "\n"
                            index = index + 1
                        except:
                            print "cannot find chains for file: " + BioUnitID.lower()
                            continue
                    out_obj.write(aline)
            except:
                continue
            # add one index for each binding site/ligand
    out_obj.close()

def aqeelNumberingParse( infile ):
    numberdict = dict()
    for line in open( infile ):
        content = line.strip().split("\t")
        numbering = content[0].strip()
        protein   = content[1].strip()
        ligand    = content[2].strip()
        proteinLigand = "\t".join( [protein, ligand] )
        numberdict[ proteinLigand ] = numbering
    return numberdict

def assignEachPDBwithNumberofMembers( PDBLeader_dict ):
    PDBMember_dict = dict()
    leader_dict = dict()
    leaderlist = []
    for each in PDBLeader_dict.keys():
        leader = PDBLeader_dict[each]
        if not leader in leaderlist:
            leaderlist.append( leader )
        try:
            leader_dict[leader] += 1
        except:
            leader_dict[leader]  = 1
    for each in PDBLeader_dict.keys():
        leader = PDBLeader_dict[ each ]
        if each in leaderlist:
            PDBMember_dict[each] = leader + "\t" + str(leader_dict[leader])
        else:
            PDBMember_dict[each] = leader + "\t" + "-"
    return PDBMember_dict

if __name__ == "__main__":
    # for everyparser of valid ligand
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )

    infiledir = "/users/ajing/ligandNet/tmp_test/final.txt"
    outfiledir1 = "/users/ajing/ligandNet/ProbisInput_test.txt"
    outfiledir2 = "/users/ajing/ligandNet/Data/IndexwithLigandInfo.txt"
    ## file for protein id : protein chains
    proteinChainFile = "/users/ajing/ligandNet/Data/proteinChain.txt"
    proteinchaindict = returnChainsForPDBID( proteinChainFile )
    ## file for aqeel's numbering for protein ligand pair
    aqeelfile = "/users/ahmedaq/work/Probis/LigandID_PDB_LigName_tab.nosql"
    aqeeldict = aqeelNumberingParse( aqeelfile )
    probisdict = convertProBiS( infiledir )
    makeProBiSInput( probisdict, everyparser.ALL, outfiledir1, outfiledir2, proteinchaindict, aqeeldict, pdb_with_numberofmembers)
