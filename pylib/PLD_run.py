#####################################
##   Purpose: for the use of Rich's binning code
##
#####################################

import sys

PLD_dir = '/users/ajing/PLD/'
sys.path.append(PLD_dir)
from getNonRedundantEntries_ajing import PLDmain

def main():
    PLDmain(90)

if __name__ == "__main__":
    main()
