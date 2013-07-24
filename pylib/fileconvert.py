"""
    fileconvert for ligand_valid.txt (output of compareMOADBiounit.py)
"""

def convertLigandValid( filedir ):
    #  convert from
    #         2lig.bio1       H_ASP   189     A
    #  to
    #          2lig.bio1      ASP.189.A
    outdir = filedir + "_new"
    outobj = open( outdir, "w" )
    for line in open(filedir):
        content = line.split("\t")
        out_line = content[0] + "\t" + content[1][2:] + "." + ".".join(content[2:-1]) + "\t" + content[-1]
        print out_line
        outobj.write( out_line )
    outobj.close()


if __name__ == "__main__":
    convertLigandValid( "/users/ajing/ligandNet/pylib/ligand_valid.txt" )
