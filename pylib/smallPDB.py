"""
    Read proteinChain.txt and produce biounit filename which contains no chain in it
    Author: ajing
    Date:   7/31/2013
"""

def getemptyFile( filedir ):
    for line in open(filedir):
        content = line.split("\t")
        print content
        if len(content) < 2:
            print line

if __name__ == "__main__":
    file = "/users/ajing/ligandNet/pylib/proteinChain.txt"
    getemptyFile(file)
