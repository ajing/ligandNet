'''
    This file is for convert original file of Rich's amino_acid_counter.pl to what ProBiS like
    Please call all files in corresponding data directory
    author: ajing
    Data  : 7/15/2013
    ChangeLog:
        7/29/2013 add another column for number of members in that family
'''

## For every_parser
import sys
from RichOutParser import RichOutParser
from Modular import PROBIS_DIR, BIOUNIT_DIR, RICHOUT_DIR, PROTEIN_DIR, MOADINDEX, OBSOLETE, EVERYCSV
__previous_pylib__ = "/users/ajing/pylib"
__INPUTDIR__ = BIOUNIT_DIR
__OUTPUT__ = PROBIS_DIR
sys.path.append(__previous_pylib__)
from every_parser import every_parser

# parameters:
alpha = 0.5   # only when ligands contains more than certain percent as the ligand name, I consider this ligand can map to certain ligand in MOAD

def processStrangeLigandName( completeLigandName ):
    # Here completeLigandName means like BGC.D or GLC BGC.B_C
    ligandName = completeLigandName.split(".")[0]
    ChainID    = completeLigandName.split(".")[1]
    ligandName.strip()
    ligandNamelist = ligandName.split()
    ChainIDlist    = ChainID.split("_")
    ligandNamelen  = len(ligandNamelist)
    if ligandNamelen == 1:
        return completeLigandName
    else:
        if len(ChainIDlist) == 1:
            ChainIDlist = ChainIDlist * ligandNamelen
        newNamelist = []
        for i in range(ligandNamelen):
            try:
                newNamelist.append( ligandNamelist[i] + "." + ChainIDlist[i] )
            except:
                print "Strange ligand name: " + completeLigandName
                raise TypeError()
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

def contains(small, big):
    if len(big) > 2:
        return contains_simple(small, big)
    if len(small) < 2:
        if small[0] in big:
            return True
        else:
            return False
    else:
        if sorted(small) == sorted(big):
            return True
        else:
            return False

def contains_simple(small, big):
    for each in small:
        if not each in big:
            return False
    if len(small) * 1.0 / len(big) > alpha:
        return True
    else:
        return False

def ligandCompare( ligand1, ligand2 ):
    list1 = sorted( ligand1.split() )
    list2 = sorted( ligand2.split() )
    list_small, list_big = sorted([list1, list2], key = len)
    if contains( list_small, list_big ):
        return True
    return False

def checkValid( PDBID, ligand, validdict ):
    ligandlist = validdict[ PDBID.upper() ]
    for each in ligandlist.keys():
        if ligandCompare( ligand, each ):
            return each
    return False

class ligandFilter:
    # only keep one ligand in output file
    def __init__(self):
        self.ligandtrack = dict()

    def checkLigand( self, PDBID, ligand ):
        if PDBID not in self.ligandtrack.keys():
            self.ligandtrack[PDBID] = [ligand]
            return False
        if ligand not in self.ligandtrack[ PDBID ]:
            self.ligandtrack[ PDBID ].append( ligand )
            return False
        return True

def readExistingOutput( infile ):
    existingcontent = []
    for line in open( infile ):
        content = line.strip().split("\t")
        aline   = "\t".join( content[1:] )
        existingcontent.append( aline )
    try:
        index = int(content[0])
    except:
        index = 1
    return ( existingcontent, index )

def indexGenerator( index ):
    while True:
        yield index
        index = index + 1

def getBindingMoadID( PDBID, ligand, MOADDict ):
    for eachligand in MOADDict[PDBID]:
        if ligandCompare( ligand, eachligand ):
            return MOADDict[PDBID][eachligand]

class oneLineInfo:
    def __init__(self):
        self.string = ""
        self.index  = 0

    def OneMoreChain( self, astring = "" ):
        # to make formate like [A: and (1,2)]
        if astring is "":
            self.string = self.string + "]"
            return
        if self.string:
            self.string = self.string + " or :" + astring
        else:
            self.string = self.string + "[:" + astring

    def addLeft( self, leftstring ):
        self.string = "\t".join( [ leftstring, self.string ] )

    def addRight( self, rightstring ):
        self.string = "\t".join( [ self.string, rightstring ] )

    def setIndex( self, index ):
        self.index = index

    def indexString( self ):
        return str( self.index ).zfill(5)

    def addIndexString( self ):
        self.string = self.indexString() + "\t" + self.string

    def writeToFile( self, fileObj ):
        # add one index for each binding site/ligand
        self.addIndexString()
        fileObj.write( self.string + "\n" )


def isleader( PDBID, PDBMemberNumberDict ):
    if PDBMemberNumberDict[PDBID][-1] == "-":
        return False
    return True

# Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
# like "[:A and (12, 13, 14, 15)]"
def makeProBiSInput( ProBiS_dict, validliganddict, outfile, outfile2, BioChainsDict = None, NumberDict = None, PDBMemberNumberDict = None ):
    # Basically in the following format:
    # Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
    # Output2 format: Index, ligandName.ligandChainID,.. ( no info for ligand residue, so ignore that )
    # Example for output2: 00001 ATP.123.A,ALA.43.C,TYR.44.C
    # 7/30/2013 only keep one ligand for one PDB
    out_obj = open( outfile, "a" )
    # 7/30/2013 Here is the structure to keep track only one case for each ligand, remove the duplication in biounit also
    onlyOneLigandEntry = ligandFilter()
    # 9/17 read existing output and generate continuous index
    existingContent, index = readExistingOutput(outfile)
    indexG = indexGenerator( index + 1 )
    for eachBioUnit in sorted( ProBiS_dict.keys() ):
        BioUnitID   = eachBioUnit
        PDBID       = BioUnitID.split('.')[0]
        # 11/12/2013 for ignoring obsolete protein
        if PDBID in OBSOLETE:
            continue
        for eachligand in sorted(ProBiS_dict[eachBioUnit].keys()):
            oneLine       = oneLineInfo()
            # The total number of binding site residue
            bindingSiteNumber = 0
            # 9/6 sort protein chain id also
            for eachproteinChainID in sorted( ProBiS_dict[eachBioUnit][eachligand].keys() ):
                ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID].sort()
                bindingSiteNumber += len(ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID])
                chainInfo = eachproteinChainID + " and (" + ",".join( map( str, ProBiS_dict[eachBioUnit][eachligand][eachproteinChainID] ) ) + ")"
                oneLine.OneMoreChain( chainInfo )
            oneLine.OneMoreChain()
            ligandName, ligandChain = eachligand.split('.')
            # 7/30/2013 for only one ligand for leader
            #if isleader(PDBID, PDBMemberNumberDict) and onlyOneLigandEntry.checkLigand( PDBID, ligandName ):
            #    continue
            if checkValid( PDBID, ligandName, validliganddict ):
                ligandChainID = processStrangeLigandName( eachligand )
                if BioChainsDict is None:
                    # whether include binding PDB chains or not
                    oneLine.addLeft("\t".join( [ BioUnitID.lower(), ligandChainID, allproteinChainInfo ] ) )
                else:
                    numbering = getBindingMoadID( PDBID.upper(), ligandName, NumberDict )
                    #print PDBID, ligandName, numbering
                    try:
                        oneLine.addLeft("\t".join( [ numbering, BioUnitID.lower(), ligandChainID, BioChainsDict[BioUnitID.lower()]]) )
                    except:
                        print [ numbering, BioUnitID.lower(), ligandChainID, BioChainsDict[BioUnitID.lower()] ]
                        raise TypeError("some keys cannot be found in a dictionary object")
                    oneLine.addRight("\t".join( [ str(bindingSiteNumber), PDBMemberNumberDict[PDBID.upper()], "".join(sorted(ProBiS_dict[eachBioUnit][eachligand].keys()))]) )
                if not oneLine.string in existingContent:
                    oneLine.setIndex( indexG.next() )
                    oneLine.writeToFile( out_obj )
    out_obj.close()

def aqeelNumberingParse( infile ):
    # extend numbering to other fields: like
    numberdict = dict()
    for line in open( infile ):
        content = line.strip().split("\t")
        numbering = content[0].strip()
        protein   = content[1].strip()
        ligand    = content[2].strip()
        if protein in numberdict:
            numberdict[ protein ][ ligand ] = numbering
        else:
            numberdict[ protein ] = dict()
            numberdict[ protein ][ ligand ] = numbering
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

def main():
    # for everyparser of valid ligand
    everyparser = every_parser(EVERYCSV)
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )
    infiledir = RICHOUT_DIR
    outfiledir2 = "/users/ajing/ligandNet/Data/IndexwithLigandInfo.txt"
    ## file for protein id : protein chains
    proteinChainFile = PROTEIN_DIR
    proteinchaindict = returnChainsForPDBID( proteinChainFile )
    ## file for aqeel's numbering for protein ligand pair
    aqeelfile   = MOADINDEX
    aqeeldict   = aqeelNumberingParse( aqeelfile )
    richout     = RichOutParser( infiledir, everyparser.ALL )
    probisdict  = richout.obj
    makeProBiSInput( probisdict, everyparser.ALL, __OUTPUT__, outfiledir2, proteinchaindict, aqeeldict, pdb_with_numberofmembers)

def test():
    # test for functions
    aqeelfile   = MOADINDEX
    aqeeldict   = aqeelNumberingParse( aqeelfile )
    print getBindingMoadID("1FPQ", "SAM", aqeeldict)
    print getBindingMoadID("1CZH", "SO4", aqeeldict)
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    ligands = "ILE LEU GLY PRO SER VAL TYR"
    print checkValid("1XR9", ligands, everyparser.ALL)
    print checkValid("4G0A", "G  U", everyparser.ALL)

if __name__ == "__main__":
    #test()
    main()
