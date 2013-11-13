"""
    This is a help python file for convertToProBiSFormat.py
    Author:ajing
    Date:  9/19/2013
"""
import sys
__previous_pylib__ = "/users/ajing/pylib"
sys.path.append(__previous_pylib__)
from every_parser import every_parser
# The original colname of file
_COLNAME = [
    "BIOUNIT",
    "BIOUNITFILE",
    "ligandName",
    "ligandChainID",
    "ligandChainIDNEW",
    "proteinChainID",
    "residueName",
    "objresidueNumber",
    "insertion",
    "AtomName",
    "AtomNumber",
    "Distance",
    "Metal"
]

class RichOutParser:
    def __init__( self, Richdir, validdict = None):
        self.infile = Richdir
        self.exceptionDIR = "/users/ajing/ligandNet/Data/joinException.txt"
        self.exceptlist = self.exceptionList()
        self.obj    = self.parseRichOutput( self.infile )
        if validdict:
            self.obj    = self.joinMultipartBindingSite( self.obj, validdict )

    def parseRichOutput( self, infile ):
        output_dict = dict()
        for line in open(infile):
            content = line.strip().split(",")
            content_dict = dict()
            tot_num = len( content )
            if tot_num != len( _COLNAME ):
                raise "Some problem with this line: " + content
            for idx in range( tot_num ):
                content_dict[ _COLNAME[idx] ] = content[ idx ]
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
                output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ].sort()
        return output_dict

    def getAllBiounitForPDB( self, PDBID, ProbisDict ):
        Biounitlist = []
        for eachBiounit in ProbisDict:
            if eachBiounit[:4] == PDBID:
                Biounitlist.append( eachBiounit )
        return Biounitlist

    def listInlist( self, smallList, largeList ):
        for each in smallList:
            if each in largeList:
                largeList.remove( each )
            else:
                return None
        return largeList

    def deepIn( self, ligand, ligandlist ):
        for each in ligandlist:
            if ligand in each:
                return each
        return False

    def updateBindingSite( self, liganddict, tmpdict ):
        for each in liganddict:
            if tmpdict.get( each, None ):
                tmpdict[ each ] = tmpdict[each] + liganddict[each]
                tmpdict[ each ] = list(set( tmpdict[ each ] ))
                tmpdict[ each ].sort()
            else:
                tmpdict[ each ] = sorted( liganddict[ each ] )
        return tmpdict

    def exceptionList( self ):
        exceptlist = []
        for eachline in open(self.exceptionDIR):
            biounit, ligandchain = eachline.strip().split("\t")
        exceptlist.append( [biounit, ligandchain] )
        return exceptlist

    def exceptionListforJoin( self, biounit, ligandandchain ):
        if [ biounit, ligandandchain ] in self.exceptlist:
            return True
        else:
            return False

    def joinforOneBiounit( self, ligands, biounit, ProbisDict):
        # return new
        # ligands is a list of few digit multipart ligand
        joinlist   = []
        ligands_cp = ligands[:]
        # get join list for ligands
        for ligandandChain in ProbisDict[biounit].keys():
            #if biounit == "2ARX.BIO1":
            #    print ligandandChain
            #    print [biounit, ligandandChain]
            if self.exceptionListforJoin( biounit, ligandandChain):
                continue
            ligand, chain = ligandandChain.split(".")
            ligandlist = ligand.split()
            ligandsleft = self.listInlist( ligandlist, ligands )
            if ligandsleft:
                joinlist.append( ligandandChain )
                ligands = ligandsleft
        if len( joinlist ) < 2 or not joinlist:
            return ProbisDict
        # reorder join list, so the same order as bindingMOAD
        ligands      = ligands_cp
        joinlist_tmp = []
        for each in ligands:
            ligand = self.deepIn( each, joinlist)
            if ligand:
                joinlist_tmp.append( ligand )
                joinlist.remove( ligand )
        joinlist = joinlist_tmp
        # join binding sites for joinlist
        # build key
        ligandlist = []
        chainlist  = []
        ProbisDict[biounit]["tmpligand"] = dict()
        for eachligand in joinlist:
            ligand, chain = eachligand.split(".")
            ligandlist = ligandlist + ligand.split()
            chainlist.append( chain )
            liganddict    = ProbisDict[biounit].pop(eachligand)
            ProbisDict[biounit]["tmpligand"] = self.updateBindingSite( liganddict, ProbisDict[biounit]["tmpligand"] )
        if ligandlist:
            key = ".".join([ " ".join(ligandlist), "_".join(chainlist) ])
            ProbisDict[biounit][key] = ProbisDict[biounit].pop( "tmpligand" )
        return ProbisDict

    def joinMultipartBindingSite( self, ProbisDict, validdict ):
        # MOADDict is a dictionary with {PDB:ligand:valid}
        for eachPDB in validdict:
            for eachligand in validdict[eachPDB]:
                ligands = eachligand.strip().split()
                if len(ligands) > 1:
                    Biolist = self.getAllBiounitForPDB( eachPDB, ProbisDict )
                    for eachBiounit in Biolist:
                        #print "eachbiounit:", eachBiounit
                        ProbisDict = self.joinforOneBiounit( ligands, eachBiounit, ProbisDict)
        #BioUnitList = self.getAllBiounitForPDB( "2ARX", ProbisDict )
        #for each in BioUnitList:
        #    print each
        #    print ProbisDict[ each ]
        return ProbisDict

if __name__ == "__main__":
    everyparser = every_parser()
    everyparser.find_PDBID_ValidLigand()
    infiledir = "/users/ajing/ligandNet/tmp_test/final.txt"
    richout = RichOutParser(infiledir, everyparser.ALL)
