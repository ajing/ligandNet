'''
    Purpose: check probis input for how many PDB and so on
    Author:  ajing
    Date:    8/21/2013
'''

__FILEDIR__ = "/users/ajing/ligandNet/pylib/ProbisInput.txt"

__COLNAME__ = [
        "BindingSiteIndex",
        "MOADIndex",
        "BioUnit",
        "LigandIdentifier",
        "ProteinChains",
        "ProbisInput",
        "BindingResidueNumber",
        "FamilyLeader",
        "FamilyMemberNumber"
        ]

def ProbisParse():
    ProbisObj = open( __FILEDIR__ )
    ProbisInfo = dict()
    # initialization
    for each in __COLNAME__:
        ProbisInfo[ each ] = []
    for line in ProbisObj:
        content = line.strip().split("\t")
        if not len(content) == len( __COLNAME__ ):
            raise "something wrong with this column: " + "\t".join( content )
        for idx in range( len(content) ):
            ProbisInfo[ __COLNAME__[idx] ].append( content[idx] )
    return ProbisInfo

def numberofPDB( Probis ):
    Biounit = Probis[ "BioUnit" ]
    PDBlist = []
    for each in Biounit:
        PDB = each.split( "." )[0]
        if not PDB in PDBlist:
            PDBlist.append( PDB )
    return len(PDBlist)

def main():
    Probisinfo = ProbisParse()
    numberPDB  = numberofPDB( Probisinfo )
    print "total number of PDBID in ProbisInput"
    print numberPDB

if __name__ == "__main__":
    main()
