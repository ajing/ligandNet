'''
    Identical ligands in biounit files
'''

from Modular import RICHOUT_DIR, RICH_COLNAME
_COLNAME = RICH_COLNAME

def getDuplicateLigandsFromRichOut( infile ):
    output_dict = dict()
    output_list = []
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
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ]  = dict()
        try:
            int(content_dict["residueNumber"])
        except:
            #print line
            #print content_dict
            continue
        if int(content_dict["residueNumber"]) in output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ]:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ] [content_dict["residueNumber"]] += 1
            outputline = ",".join([content_dict["BIOUNIT"], content_dict["ligandName"], content_dict["ligandChainID"], content_dict["residueNumber"]])
            print outputline
            output_list.append(outputline)
        else:
            output_dict[ content_dict["BIOUNIT"] ][ content_dict["ligandName"] + "." + content_dict["ligandChainID"] ][ content_dict["proteinChainID"] ] [content_dict["residueNumber"]] = 1



if __name__ == "__main__":
    getDuplicateLigandsFromRichOut(RICHOUT_DIR)
