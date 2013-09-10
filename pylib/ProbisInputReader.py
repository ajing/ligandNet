"""
    This file is to read probis input and do some statistics
    Author: ajing
    Date:   9/10/2013
"""

class Probis:
    def __init__( self ):
        self.PROBISFILE = "/users/ajing/ligandNet/pylib/ProbisInput.txt"
        self.COLNAME = [
                "INDEX",
                "MOADINDEX",
                "BIOUNIT",
                "LIGANDID",
                "PROTEINCHAIN",
                "BINDINGSITE",
                "BINDINGSITESIZE",
                "HEADER",
                "ISHEADER"
                ]
        self.probisdict = self.ProbisInputParser( self.PROBISFILE )

    def ProbisInputParser( self, infile ):
        ProbisDict = dict()
        for each in self.COLNAME:
            ProbisDict[ each ] = []
        for line in open(infile):
            content = line.strip().split("\t")
            total_len = len(content)
            if not total_len == len( self.COLNAME):
                raise "file corrupted as different colnames"
            for each in range( total_len ):
                ProbisDict[ self.COLNAME[each] ].append( content[each] )
        return ProbisDict

if __name__ == "__main__":
    test = Probis()
