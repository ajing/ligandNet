#################################
##  Author: ajing
##  Date:   6/20/2013
##  Purpose:For rich's code only take four letters as PDB indetification number, so I map a whole directory to a four letter "*.pdb"
##  I/O:    input is a directory and output is a new directory with renamed files and also a hashtable
#################################

import shutil
import os

def renameFile( filename, number, olddir, newdir ):
    src = os.path.join( olddir, filename )
    newfilename = hex(number).split("x")[-1].zfill(4) + ".pdb"
    dst = os.path.join( newdir, newfilename )
    shutil.copyfile( src, dst )
    return ( filename, newfilename )

def goThroughAllfiles( sourcedir, destdir, hashfile ):
    index = 1
    hashfile_obj = open( hashfile, "w" )
    for each in os.listdir(sourcedir):
        print each
        source, destination = renameFile( each, index, sourcedir, destdir )
        hashfile_obj.write( source + "\t" + destination + "\n" )
        index += 1

def main():
    srcdir = "/users/ajing/BindingMOAD/biou_unzip"
    dstdir = "/users/ajing/BindingMOAD/biou_rename"
    hfile  = "/users/ajing/BindingMOAD/hashtable.txt"

    if not os.path.exists( srcdir ):
        print "cannot find source directory"
    if not os.path.exists( dstdir ):
        os.makedirs( dstdir )
    goThroughAllfiles( srcdir, dstdir, hfile )

if __name__ == "__main__":
    main()
