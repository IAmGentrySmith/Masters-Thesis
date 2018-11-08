from chimera import runCommand as rc
from chimera import replyobj
import sys
import os

#standard sys.argv[] for script args?
# sys.argv[0] = directory
os.chdir(sys.argv[0])

file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]
fn = file_names[0]
# inPDB = chimera.openModels.open('/Users/gentry/Desktop/test/testmol.pdb', type="PDB")

rc("open " + fn)

rc("rotation 1 reverse #0:1.HET@/serialNumber=2 #0:1.HET@/serialNumber=3")

for i in range(72):
    #replyobj.status("Processing " + fn)
    #rc("open " + fn)
    #rc("rotation 1 reverse #0:1.HET@/serialNumber=2 #0:1.HET@/serialNumber=3")
    rc("rotation 1 5")
    newName = (fn[:-3] + str((i*5)) + ".pdb")
    rc("write format pdb 0 " + newName)
    #rc("close ")


    # chimera.runCommand("rotation 2 3 5")
    # newName = ( inPDB[:-3] + i*5 + ".pdb" )
    # chimera.runCommand("write format pdb " + newName)

