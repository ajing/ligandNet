'''
    Purpose: check proteinChain.txt to remove some empty entries
    Author: ajing
    Date:   8/21/2013
'''

__PROTEINCHAIN__ = "/users/ajing/ligandNet/pylib/proteinChain.txt"

def removeEmpty( filedir ):
    f = open( filedir, 'r+')
    newText = []
    noChainIDnum = 0
    for line in f:
        content = line.strip().split("\t")
        print content
        if len(content) > 1:
            newText.append(line)
        else:
            noChainIDnum = noChainIDnum + 1
    print "Total number PDB has no Chain ID: " + str( noChainIDnum )
    f.seek(0)
    f.write("".join(newText))
    f.truncate()
    f.close()

def main():
    removeEmpty( __PROTEINCHAIN__ )

if __name__ == "__main__":
    main()
