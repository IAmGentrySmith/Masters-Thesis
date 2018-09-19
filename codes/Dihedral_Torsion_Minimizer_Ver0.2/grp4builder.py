#!/usr/bin/python

import sys
import subprocess

# argument: sys.argv[num]
# Replacement: sed -i -e 's/IN/OUT/g' FILE > NEWFILE

inFile = file(sys.argv[1])

def DoIT():
    for first in {' C', 'Si', 'Ge'}:
        name1 = "%s" % (first.lstrip(' '))
        out1 = open(name1, "w")
        cmdStr = "sed -e 's/1 GE/1 %s/g' ./%s >> ./%s.pdb" % (first, inFile, name1)
        # subprocess.call(cmdStr, shell=True, stdout=out1)
        subprocess.Popen(cmdStr, shell=True, executable='/bin/bash')
        out1.close()
        for second in {' C', 'Si', 'Ge'}:
            name2 = name1 + "_%s" % (second.lstrip(' '))
            out2 = open(name2, "w")
            cmdStr = "sed -e 's/2 GE/2 %s/g' ./%s.pdb >> ./%s.pdb" % (second, name1, name2)
            # subprocess.call(cmdStr, shell=True, stdout=out2)
            subprocess.Popen(cmdStr, shell=True, executable='/bin/bash')
            out2.close()
            for third in {' C', 'Si', 'Ge'}:
                name3 = name2 + "_%s" % (third.lstrip(' '))
                out3 = open(name3, "w")
                cmdStr = "sed -e 's/3 GE/3 %s/g' ./%s.pdb >> ./%s.pdb" % (third, name2, name3)
                # subprocess.call(cmdStr, shell=True, stdout=out3)
                subprocess.Popen(cmdStr, shell=True, executable='/bin/bash')
                out3.close()
                for fourth in {' C', 'Si', 'Ge'}:
                    name4 = name3 + "_%s" % (fourth.lstrip(' '))
                    out4 = open(name4, "w")
                    cmdStr = "sed -e 's/4 GE/4 %s/g' ./%s.pdb >> ./%s.pdb" % (fourth, name3, name4)
                    # subprocess.call(cmdStr, shell=True, stdout=out4)
                    subprocess.Popen(cmdStr, shell=True, executable='/bin/bash')
                    out4.close()

DoIT()