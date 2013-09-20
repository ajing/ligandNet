'''
    This file is for convert original file of Rich's amino_acid_counter.pl to what ProBiS like
    author: ajing
    Data  : 7/15/2013
    ChangeLog:
        7/29/2013 add another column for number of members in that family
'''

## For every_parser
import sys
import os
from RichOutParser import RichOutParser
__previous_pylib__ = "/users/ajing/pylib"
#__INPUTDIR__ = "/users/ajing/ligandNet/2012_biounits/"
__INPUTDIR__ = "/users/ajing/ligandNet/BindingMoad2011_test/BindingMoad2011/"
__OUTPUT__ = "/users/ajing/ligandNet/ProbisInput_test.txt"
sys.path.append(__previous_pylib__)
from every_parser import every_parser

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
                raise TypeError
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
    for i in xrange(len(big)-len(small)+1):
        for j in xrange(len(small)):
            if big[i+j] != small[j]:
                break
        else:
            return i, i+len(small)
    return False

def ligandCompare( ligand1, ligand2 ):
    list1 = sorted( ligand1.split() )
    list2 = sorted( ligand2.split() )
    if contains( list1, list2 ) or contains( list2, list1 ):
        return True
    return False

def checkValid( PDBID, ligand, validdict ):
    ligandlist = validdict[ PDBID.upper() ]
    for each in ligandlist.keys():
        if ligandCompare( ligand, each ):
            return ligandlist[ each ]
    return False

class ligandFilter:
    # only keep one ligand in output file
    def __init__(self):
        self.ligandtrack = dict()

    def checkLigand( self, PDBID, ligand ):
        if PDBID not in self.ligandtrack.keys():
            self.ligandtrack[PDBID] = []
            return False
        if ligand not in self.ligandtrack[ PDBID ]:
            self.ligandtrack[ PDBID ].append( ligand )
            return False
        return True

def readExistingOutput():
    existingcontent = []
    for line in open( __OUTPUT__ ):
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

    def indexAdd( self ):
        self.index = indexGenerator()

    def indexString( self ):
        return str( self.index ).zfill(5)

    def addIndexString( self ):
        self.string = self.indexString() + "\t" + self.string

    def writeToFile( self, fileObj ):
        self.addIndexString()
        fileObj.write( self.string + "\n" )


# Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
# like "[:A and (12, 13, 14, 15)]"
def makeProBiSInput( ProBiS_dict, validliganddict, outfile, outfile2, BioChainsDict = None, NumberDict = None, PDBMemberNumberDict = None ):
    # Basically in the following format:
    # Output format: Index, BioUnitID, ligandChainID, [:proteinChainID and residueNumber]
    # Output2 format: Index, ligandName.ligandChainID,.. ( no info for ligand residue, so ignore that )
    # Example for output2: 00001 ATP.123.A,ALA.43.C,TYR.44.C
    # 7/30/2013 only keep one ligand for one PDB
    out_obj = open( outfile, "a" )
    ExistingPairs = []
    # 7/30/2013 Here is the structure to keep track only one case for each ligand, remove the duplication in biounit also
    onlyOneLigandEntry = ligandFilter()
    # 9/17 read existing output and generate continuous index
    existingContent, index = readExistingOutput()
    indexG = indexGenerator( index )
    for eachBioUnit in ProBiS_dict.keys():
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
            ligandName  = eachligand.split('.')[0]
            ligandChain = eachligand.split('.')[1]
            # 7/30/2013 for only one ligand
            #if onlyOneLigandEntry.checkLigand( PDBID, ligandName ):
            #    continue
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

def existingStatistics( probisdict ):
    biounits = probisdict.keys()
    print "Number of biounit files: " + str(len(biounits))
    PDBs     = set( [ each.split(".")[0] for each in biounits ])
    print "Number of PDBs: " + str(len(PDBs))
    filelist = os.listdir( __INPUTDIR__ )
    print "Total number of biounit files: " + str(len(filelist))
    AllPDBs     = set( [ each.split(".")[0] for each in filelist])
    print "Total number of PDBs: " + str(len(AllPDBs))

def main():
    # for everyparser of valid ligand
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    pdb_with_numberofmembers = assignEachPDBwithNumberofMembers( everyparser.ALL_leader )

    infiledir = "/users/ajing/ligandNet/tmp_test/final.txt"
    outfiledir2 = "/users/ajing/ligandNet/Data/IndexwithLigandInfo.txt"
    ## file for protein id : protein chains
    proteinChainFile = "/users/ajing/ligandNet/Data/proteinChain.txt"
    proteinchaindict = returnChainsForPDBID( proteinChainFile )
    ## file for aqeel's numbering for protein ligand pair
    aqeelfile = "/users/ahmedaq/work/Probis/LigandID_PDB_LigName_tab.nosql"
    aqeeldict   = aqeelNumberingParse( aqeelfile )
    richout     = RichOutParser( infiledir, everyparser.ALL )
    probisdict  = richout.obj
    makeProBiSInput( probisdict, everyparser.ALL, __OUTPUT__, outfiledir2, proteinchaindict, aqeeldict, pdb_with_numberofmembers)

def main2():
    infiledir = "/users/ajing/ligandNet/tmp_test/final.txt"
    probisdict = parseRichOutput( infiledir )
    existingStatistics( probisdict )

if __name__ == "__main__":
    main()
    #main2()

