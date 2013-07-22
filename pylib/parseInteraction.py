###############################
##   Date  : 06/18/2013
##   Author: ajing
##   Purpose: parse ligand interaction output from Rich
##   I/O   : input is the output of amino_acid_counter.pl
###############################

from collections import namedtuple

ligandInteraction = namedtuple( "ligandInter", "ligandname, chainid")

def changeFileFormate( filedir, tmpdir ):
    # add first column as PDBID and strip some space
    PDBID = filedir.split("/")[-1][:4]
    print PDBID
    tmpobj = open( tmpdir, "w")
    for line in open(filedir):
        content = line.strip().split("\t")
        newcontent = []
        for each in content:
            newcontent.append(each.strip())
        newcontent.pop(5)
        newcontent.insert( 0, PDBID.upper() )
        tmpobj.write( ",".join(newcontent) + "\n")

if __name__ == "__main__":
    fileDIR = "/users/ajing/perllib/1hg0.out"
    changeFileFormate( fileDIR, "tmp.txt")

