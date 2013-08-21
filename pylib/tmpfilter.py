filedir = "/users/ajing/ligandNet/pylib/emptyFileName"
outdir  = filedir + "_convert"

PDBlist = []

out_obj = open(outdir, "w")
for line in open(filedir):
    PDB = line.strip().split(".")[0]
    if not PDB in PDBlist:
        out_obj.write(line)
        PDBlist.append("")
