from Bio.Blast import NCBIStandalone

blastdir = "/users/ajing/BindingMOAD/blastout/1gkc_a.out"

def printAllAttributes(object):
    for each in dir(object):
        print each
        if each[:2] != "__":
            print getattr( object, each )
    # From here I know query and sbjct is what I want

def printAtt( blast_record ):
    for alignment in blast_record.alignments:
        #print dir(alignment)
        print alignment.title
        print alignment.length
        for hsp in alignment.hsps:
            print hsp
            printAllAttributes(hsp)
            break

#def mapQuery( hsp_query ):

if __name__ == "__main__":
    result_handle = open( blastdir )
    blast_parser = NCBIStandalone.BlastParser()
    blast_record = blast_parser.parse(result_handle)
    print dir(blast_record)
    printAtt(blast_record)

