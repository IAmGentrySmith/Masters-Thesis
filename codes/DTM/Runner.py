### Author: Gentry Smith
### Date: April 22, 2017
### Description: This is the runner file that is the primary executable for the torsion minimizer. Currently is the
###     only file utilized.

# Inputs:
# Arg1: the molecule file to be minimized (currently only accepts a pdb file)

import sys
import subprocess
import math

# IO Validator: validates user-submitted molecule.
def IOValidator():
    isValid = False
    # Check for valid length of args (2)
    if len(sys.argv) == 2:
        # Check arg to make sure it's a file.
        argFile = sys.argv[1]
        try:
            inputFile = open(argFile)
            # Finally, make sure the file is a .pdb
            if inputFile[-4:] == ".pdb":
                isValid = True
            else:
                print("This is not a .pdb file. Please try again with a .pdb file.\n")
            inputFile.close()
        except IOError:
            print("System was not able to open '", str(argFile), "'.")
    # too long
    elif len(sys.argv) > 2:
        print("You have too many arguments. Call the file as 'Runner.py [molecule file]' and try again.\n")
    # too short
    else:
        print("You do not have enough arguments. Start the program as 'Runner.py [molecule file]' and try again.\n")
    # return validity boolean
    return isValid

# Get Torsions: initiates function to get user-specified torsion bonds. Returns bonds as int[[a,b],[a,b]] list
def getTorsions():
    torsions = [[0, 0]]
    newTorsion = "first"
    firstTime = True
    doneCheck = ""
    badIn = False

    # loop for all torsions until user types "done"
    while newTorsion != "":
        if firstTime:
            print("It's time to define the torsions of the molecule and declare which bonds you would like to rotate.\n")
            print("Before going any further, it's important to note at this time that version 0.2 (current) will assume the torsions you enter are completely correct. You'll see a bunch of error messages soon if it isn't correct.\n")
            print("Open the .pdb file and identify the numbers of the atoms on the .pdb that will make the bond (the first number on the line of each atom)\n\n")
            print("Now it's time to enter in the numbers of the two atoms. We'll do it one at a time.")

            firstTor = raw_input("Type in the number of the first atom in the bond and hit enter. \nEx: type 3 and then hit enter.\n")

            try:
                confFirstTor = int(firstTor)
            except ValueError:
                print("You typed in '", firstTor, "', which is not a number. Let's start again.")
                badIn = True

            secondTor = raw_input("Type in the number of the second atom in the bond and hit enter. \nEx: type 3 and then hit enter.\n")

            try:
                confSecondTor = int(secondTor)
            except ValueError:
                print("You typed in '", secondTor, "', which is not a number. Let's start again.")
                badIn = True
            firstTime = False

        else:
            print("Open the .pdb file and identify the numbers of the atoms on the .pdb that will make the bond (the first number on the line of each atom)\n\n")

            firstTor = raw_input("Type in the number of the first atom in the bond and hit enter. \nEx: type 3 and then hit enter.\n")

            try:
                confFirstTor = int(firstTor)
            except ValueError:
                print("You typed in '", firstTor, "', which is not a number. Let's start again.")
                badIn = True

            secondTor = raw_input("Type in the number of the second atom in the bond and hit enter. \nEx: type 3 and then hit enter.\n")

            try:
                confSecondTor = int(secondTor)
            except ValueError:
                print("You typed in '", secondTor, "', which is not a number. Let's start again.")
                badIn = True
            firstTime = False

        if badIn == False:
            newTorsion  = [confFirstTor, confSecondTor]
            if torsions == [[0, 0]]:
                print("You added a new torsion: ", newTorsion, "\n")
                torsions = newTorsion
            else:
                torsions.append(newTorsion)
                print("The current torsions you have created are:\n")
                for each in torsions:
                    print(each, "\n")
            doneCheck = raw_input("If you would like to add another torsion, press enter. If you are finished adding torsions, type 'done' and press enter\n")

            if str(doneCheck) == "done":
                print("Finished entering torsions. Begining the work.\n")
            else:
                newTorsion = "first"

        if badIn == True:
            firstTime = True
            badIn = False
            newTorsion = "first"

    return torsions

# Get Conformation Count: determines conformations needed. Returns list in form: [#conf, rotDeg, rotRng]
def getConformationInfo(depth, torsions):
    # rotates 60 degrees on the first search, then logarithmic decrease from 10 for each subsequent search.
    rotDeg = [60, 10]
    # full torsion range for first search, logarithmic decrease from 50 for each subsequent search
    rotRng = [360, 50]
    # number of conformations needed
    numConf = 0
    # degrees per rotation
    deg = 0
    # rotation range
    rng = 0
    # number of rotations per torsion
    rotTick = 0

    # determine counts from depth
    if depth >= 2:
        deg = math.pow(10, (2-depth))
        rng = deg*5
    elif depth <2:
        deg = rotDeg[depth]
        rng = rotRng[depth]
    if depth == 1:
        rotTick = 6
    elif depth >= 1:
        rotTick = 11

    numConf = math.pow(torsions, rotTick)

    return [numConf, deg, rng]



def Launcher():
    valid = IOValidator()
    if valid:
        # do everything
        depth = 0

        InitWD()


    else:
        print("There was a problem while reading in the molecule file. Please try again.\n")
        exit()


# Initiates proper working directory.
def InitWD():


# Recursive search through molecule torsions
def RecursiveSearch(depth):

    torsions = getTorsions()


Launcher()