#!/usr/bin/env pychimera

# Author: Gentry Smith, Oklahoma State University
# Date last modified: June 20, 2017

### This standalone file utilizes PyChimera to enable Python-managed commands to utilize Chimera via CLI.
###   Currently, this file only evaluates the desired torsions and builds them all to a directory
###   created with the same name as the passed molecule.

# args passed: ChimeraRunner.py [rotCmd] [rotDist] [rotTicks] [torsionData] [inputPDB]
# rotCmd: full, both (full:all the way around. both: rotate dist/ticks in both directions {i.e. expanding 0 +/- 2}
# rotDist:
# rotTicks:
# torsionData:
# inputPDB:

import sys
import subprocess

from chimera import *
my_mod=chimera.openModels.open('1crn', type="PDB")