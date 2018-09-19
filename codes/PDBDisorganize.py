
#!/usr/bin/python

# Author = Gentry Smith
# Copyright 2016, all rights reserved

# this reads in a .PDB file, takes an argument for deformities per molecules, and randomly organizes the crystal
# structure into a disordered proton formation

# import sample: python PDBDisorganize.py arg1 arg2 arg3
# where:
# arg1 = source pdb file to be read (ex: acetone.pdb or acetone)
# arg2 = number of defects per molecule (in H20, num of non-hydrogen-bonds. from 0 to 4)
# arg3 = desired output pdb file name

import sys
print sys.path
import string
import numpy as np
import math
import random

sys.setrecursionlimit(10000000) # maximum recursive depth. Set to (10,000,000) as under maximum


pdbIN = file(sys.argv[1])     # source PDB file
maxErr = int(sys.argv[2])    # max errors allowed
pdbOUT = str(sys.argv[3])    # output file name
finalData = [ [ [ 0 for i in range(3) ] for j in range(3) ] for k in range(300) ]

# looks at args validity
def checkArgs(arg1, arg2, arg3):
    returnBool = False
    if type(arg1) != file: # check arg1
        print"Bad arg", arg1, " must be a file "
        returnBool = True
    if type(arg3) != str: # check arg3
        print"Bad arg", arg3, ", must be a file name"
        checkPDBSuffix(arg3)
        print arg3
        returnBool = True
    if type(arg2) != int:  # check arg2 type
        print "Bad arg2: ", arg2, " is not an int."
        returnBool = True
    elif type(arg2) == int:
        if arg2 < 0 or arg2 > 4:  # check arg2 range
            print "arg2 is not in a valid range 0 <= arg2 <= 4"
            returnBool = True
    return returnBool

def checkPDBSuffix(pdbFile):
    if string.find(pdbFile, '.pdb', 0, len(pdbFile)) == -1:
        print("did not find 'pdb' in ", pdbFile, ". Appending...")
        pdbFile += '.pdb'



# reads in file,
def readFile(fileName):
    print "Reading file..."
    # gets number of atoms
    atoms = 0
    for line in fileName:
        data = line.split()
        if len(data) > 0:
            if data[0] != "CONECT" and data[0] != "END":
                atoms += 1
    # print "atoms: ", atoms
    numMol = atoms / 3  # assumes 3-atom water molecule
    dataTable = [ [ [ 0 for i in range(3) ] for j in range(3) ] for k in range(numMol) ]
    fileName.seek(0)
    iter0 = 0
    iter1 = 0
    pdbType = -1
    for line in fileName:
        data = line.split()
        if pdbType == -1:
            if data[0] == "ATOM":
                pdbType = 0
            elif data[0] == "HETATM":
                pdbType = 1
        # print "LineTuple= ", data
        if len(data) > 1 and ( data[0] == "ATOM" or data[0] == "HETATM" ):
            if data[0] == "ATOM":
                newData = getDataATOM(data)
                for i in range(3):
                    #data[molecule][atom][X/Y/Z]
                    dataTable[iter0][iter1 % 3][i] = newData[i]
            elif data[0] == "HETATM":
                dataTable[iter0][iter1 % 3] = getDataHETATM(data)
            if iter1 == 2:
                iter0 += 1
                iter1 = 0
            elif iter1 != 2:
                iter1 += 1
    # print "DataTable: ", dataTable
    print "File read"
    return dataTable, pdbType


    # Split by index
    # if having a problem with reading data, check .pdb to see if data has a space between each value

# reads XYZ coordinate data from ATOM-type pdb
def getDataATOM(strLine):
    # print "Getting ATOM Data..."
    dataLine = strLine[5:8]
    # print "dataline: ", dataLine
    i = 0
    while i < 3:
        # print "dataline[", i, "]: ", dataLine[i]
        dataLine[i] = float(dataLine[i])
        # print "dataline[", i, "] type: ", type(dataLine[i])
        i += 1
    return dataLine


# reads XYZ coordinate data from HETATM-type pdb
def getDataHETATM(strLine):
    # print "Getting HETATM Data..."
    dataLine = strLine[5:8]
    # print "dataline: ", dataLine
    i = 0
    while i < 3:
        # print "dataline[", i, "]: ", dataLine[i]
        dataLine[i] = float(dataLine[i])
        # print "dataline[", i, "] type: ", type(dataLine[i])
        i += 1
    return dataLine


# gets all four position vectors of hydrogen/lone pair as offset of oxygen molecule
def getOrientations( molecule ):
    # 120 degrees = ( 2 * pi ) / 3 radians
    theta = ( ( 2 * math.pi ) / 3 )
    newMol = zeroOrientation(molecule)
    returnInt1 = rotateMolecule(newMol[1], newMol[2], theta)
    returnInt2 = rotateMolecule(newMol[1], newMol[2], (-1 * theta) )
    return [returnInt1, returnInt2]


# randomly selects new orientation, returns two unique ints, from 0 to 3 inclusively
def newRandOrientation( positions ):
    # print "Changing orientation"
    randVal1 = random.randint(0,3)
    randVal2 = random.randint(0,3)
    while randVal1 == randVal2:
        randVal2 = random.randint(0,3)
    newMol = [ [ 0, 0, 0 ],
                positions[ randVal1 ] ,
                positions[ randVal2 ]  ]
    return newMol

# selects new orientation from list. Reduces computational overhead in re-orientation option traversal
def newSetOrientation( positions, pos1, pos2 ):
    newMol = [ [ 0, 0, 0 ],
               positions[ pos1 ],
               positions[ pos2 ] ]
    return newMol


# sets molecule coordinates so that oxygen is the origin
def zeroOrientation(source):
    # print "Zeroing Molecule..."

    oxy = source[0]
    hyd1 = source[1]
    hyd2 = source[2]

    # print "Oxygen pos: ", oxy
    # print "Hydrogen 1: ", hyd1
    # print "Hydrogen 2: ", hyd2

    zeroedOrigin = [0, 0, 0]
    zeroedHyd1 = [0, 0, 0]
    zeroedHyd2 = [0, 0, 0]
    for i in range(3):
        zeroedHyd1[i] = hyd1[i] - oxy[i]
        zeroedHyd2[i] = hyd2[i] - oxy[i]

    # print "Zeroed Hydrogen 1: ", zeroedHyd1
    # print "Zeroed Hydrogen 2: ", zeroedHyd2

    # return new molecule position
    newMol = [zeroedOrigin, zeroedHyd1, zeroedHyd2]
    return newMol

# resets the zeroed molecule to the original oxygen position
def resetOrientation(oxygenPos, molecule):
    # print "Resetting molecule..."
    rO = oxygenPos
    rH1 =[0,0,0]
    rH2 =[0,0,0]
    newMol = []
    for i in range(3):
        rH1[i] = molecule[1][i] + rO[i]
        rH2[i] = molecule[2][i] + rO[i]
        newMol = [rO, rH1, rH2]
    # print "Rebuilt Molecule: ", newMol
    return newMol

# rotates vector about axis for theta degrees
# Handler for rotationMatrix function below
def rotateMolecule(vector, axis, theta):
    rotMatx = rotationMatrix(axis, theta)
    return np.dot(rotMatx, vector)


# Creates Rotation matrix for a given axis and theta
# from stackoverflow user unutbu
# page: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
def rotationMatrix(axis, theta):
    """

    :type axis: list
    :type theta: union
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis /= math.sqrt(np.dot(axis, axis))
    a = math.cos( theta/2.0 )
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = (a * a), (b * b), (c * c), (d * d)
    bc, ad, ac, ab, bd, cd = (b * c), (a * d), (a * c), (a * b), (b * d), (c * d)
    return np.array( [ [ (aa + bb - cc - dd), ( 2 * ( bc + ad ) ), ( 2 * ( bd - ac ) ) ],
                       [ ( 2 * ( bc - ad ) ), (aa + cc - bb - dd), ( 2 * ( cd + ab ) ) ],
                       [ ( 2 * ( bd + ac ) ), ( 2 * ( cd - ab ) ), (aa + dd - bb - cc) ] ] )


# gets results from rotateAboutAxis plus two Hydrogens to get the tetrahedron positions
def getTetrahedronPositions(molecule):
    positions = [ [ 0 for i in range(3) ] for j in range(4) ]
    newMol = zeroOrientation(molecule)  # zero molecule
    positions[0] = newMol[1]
    positions[1] = newMol[2]
    newPos = getOrientations(molecule)  # get final two positions
    positions[2] = list(newPos[0])
    positions[3] = list(newPos[1])
    return positions                    # return all four positions


# checks distance of new positions from zero
def checkDist(posArray):
    distance = [0 for i in range(len(posArray))]
    for i in range(len(posArray)):
        distance[i] =  ( (posArray[i][0] * posArray[i][0]) +
                         (posArray[i][1] * posArray[i][1]) +
                         (posArray[i][2] * posArray[i][2]) )
        # print "Distance", i, ": ", distance[i]
    avg = 0
    for i in range(len(posArray)):
        avg += distance[i]
    averageDistance = ( avg / len(posArray) )
    # print "Average Distance: ", averageDistance
    return averageDistance


# prints data given a 3D table of water molecules
def printData(data):
    print "Data: "
    strData = [" O", "H1", "H2"]
    dimData = ["X", "Y", "Z"]
    bigAvg = 0
    numAtoms = 0
    for mol in range(len(data)):
        for atom in range(len(data[mol])):
            printStr = str(mol) + ": " + strData[atom] + ": "
            for dimension in range(3):
                printStr += dimData[atom] + ":" + "{:7.3f}".format(data[mol][atom][dimension]) + "\t"
            print printStr
        bigAvg += checkDist(zeroOrientation(data[mol])[1:])
        numAtoms += 1
        print ""
    print "total average distance: ", bigAvg / numAtoms


# checks validity of molecule
def isDefectiveCheck(err, neighborData, posData, index):
    # find nearby molecules (avg oxygen distance???)
    print "checking for defects at index", index, "..."
    print "neighbor indices: ", neighborData[index]
    returnBool = False
    neighbors = 4
    for i in range(4):  # count real neighbors
        if neighborData[index][1][i] == -1:
            neighbors -= 1
    if neighbors <= err:    # de facto good if num(neighbors) < maxErrAllowed
        # print "Fewer neighbors than allowed errors. de facto Good Orientation"
        returnBool = True
    elif neighbors > err:   # enough neighbors to require check
        # print "More neighbors than error threshold"
        defectCount = 0
        for neighbor in range(4):   # check each neighbor
            if neighborData[index][1][neighbor] != -1:  # skip over non-existant neighbors
                molA = posData[index]
                molB = posData[ neighborData[index][1][neighbor] ]
                oxyDist = getDistBetweenAtoms(molA[0], molB[0])

                if minHydrogenDistance(molA, molB) > oxyDist: # check for facing lone pairs
                    print "Double Lone Pair defect"
                    defectCount += 1
                    break
                else:   # check for facing protons
                    smallerHydrogenDistanceCount = 0
                    isDefective = False
                    for first in range(2):
                        if not isDefective:
                            for second in range(2):
                                newDist = getDistBetweenAtoms(molA[first + 1], molB[second + 1])
                                if newDist < oxyDist:
                                    smallerHydrogenDistanceCount += 1
                            if smallerHydrogenDistanceCount > 1:
                                print "Double Hydrogen defect"
                                defectCount += 1
                                isDefective = True
        # print "Defects found:", defectCount
        if defectCount > 4:
            print "IMPOSSIBLE AMOUNT OF DEFECTS DETECTED!!!! AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
        if defectCount > err:
            # print "Found a bad molecule!"
            returnBool = False
        else:
            # print "Molecule is within parameters."
            returnBool = True

    return returnBool


# randomly re-reorients molecule and neighbors, rechecks all
def rerunMolAndNeighbors(err, neighborData, posData, index):
    # print "Re-reordering molecule at", index
    # err - max errors allowed
    # neighborData - int[4] of neighbor indices
    # posData - array of all molecule position vectors
    # index - location of focus molecule in posData
    isGood = False
    timeCount = 0
    while not isGood:
        # re-rotate molecule through all positions (iterated through all orientations)
        positions = getTetrahedronPositions(posData[index])
        zeroedMol = newRandOrientation(positions)
        # print "isGood CHECK", isGood
        isGood, posData = iterThroughRotations(err, neighborData, posData, index)
        posData[index] = resetOrientation(posData[index][0], zeroedMol)
        if timeCount >= 13: # { (1 - 1/6)^n < 0.05 } says n = 17
            # BROKEN - need to rebuild
            # 0. evaluated molecule has too many defects
            # 1. reorient molecule statistically probable amount of times to cover all orientations
            # 2. Repeat 1. with neighbor 1
            # 2a repeat 1. with original molecule
            # 3. Repeat 2. with neighbor 2, 3, 4, as/if necessary
            for neighborIndex in range(4):
                if neighborData[index][1][neighborIndex] != -1:
                    positions = getTetrahedronPositions(posData[neighborIndex])
                    zeroedMol = newRandOrientation(positions)
                    posData[neighborIndex] = resetOrientation(posData[neighborIndex][0], zeroedMol)
                    # isGood = isDefectiveCheck(err, neighborData, posData, neighborIndex)
            isGood = isDefectiveCheck(err, neighborData, posData, index)
            if not isGood:
                isGood, posData = rerunMolAndNeighbors(err, neighborData, posData, neighborData[index][1][neighborIndex])
    finalData = posData
    return True, finalData

# iterates molecule through all possible rotations
def iterThroughRotations(err, neighborData, posData, index):
    isGood = False
    pos1 = 0  # tetrahedral position for H1
    pos2 = 0  # tetrahedral position for H2
    while not isGood or (pos1 != 3 and pos2 != 3):  # iterates through all orientations, stops if good orientation
        if pos1 != pos2:
            posData[index] = newSetOrientation(posData[index][0], pos1, pos2)
            isGood = isDefectiveCheck(err, neighborData, posData, index)
        if pos2 < 3:
            pos2 += 1
        elif pos2 == 3:
            if pos1 < 3:
                pos1 += 1
                pos2 = 0
    return isGood, posData
# determines minimum hydrogen distance between two atoms
def minHydrogenDistance(mol1, mol2):
    minDist = 100
    for first in range(2):
        for second in range(2):
            newDist = getDistBetweenAtoms(mol1[first+1], mol2[second+1])
            if newDist < minDist:
                minDist = newDist
    return minDist





# finds neighboring molecules of each molecule
def getNeighbors(data):
    returnData = [ [ [ 0 for i in range(4) ] for j in range(2) ] for k in range(len(data)) ] # data[molecule][distance,index][four values]
    for mol1 in range(len(data)):
        minDist = [100, 100, 100, 100]
        minIndex = [0, 0, 0, 0]
        for mol2 in range(len(data)):
            if mol1 != mol2:
                newMin = getDistBetweenAtoms(data[mol1][0], data[mol2][0])

                bigIndex = indexOfBiggest(minDist)
                if newMin < minDist[bigIndex]:
                    minDist[bigIndex] = newMin
                    minIndex[bigIndex] = mol2
        for i in range(4):
            if minDist[i] >= 9:
                minDist[i] = -1
                minIndex[i] = -1
        # print "Four smallest Distances of", mol1, ": ", minDist
        # print "Four smallest Indices of", mol1, ": ", minIndex
        returnData[mol1] = [minDist, minIndex]
    return returnData


# finds distance between oxygen atoms
def getDistBetweenAtoms( mol1, mol2 ):
    distance = ( ( ( mol1[0] - mol2[0] ) * ( mol1[0] - mol2[0] ) ) +
                 ( ( mol1[1] - mol2[1] ) * ( mol1[1] - mol2[1] ) ) +
                 ( ( mol1[2] - mol2[2] ) * ( mol1[2] - mol2[2] ) ) )
    return distance

# gets index of largest item from a list
def indexOfBiggest(check):
    bigIndex = 0
    for i in range(len(check)):
        if check[i] > check[bigIndex]:
            bigIndex = i
    return bigIndex


# writes data to PDB file
def writeDataPDB(data, pdbType):
    print "Writing Data to", str(pdbOUT)
    fileName = str(pdbOUT)
    output = open(fileName, 'w')
    if pdbType == 0:
        writeDataPDBATOM(data, output)
    elif pdbType == 1:
        writeDataPDBHETATM(data, output)
    output.close()


# Writes data to PDB file style = ATOM
def writeDataPDBATOM(data, inFile):
    iterator = 0
    for molecule in range(len(data)):
        for atom in range(3):
            iterator += 1
            outStr = "ATOM  "
            outStr += str(iterator)
            while len(outStr) < 11:
                outStr = outStr[:6] + " " + outStr[6:]
            outStr += " "
            if atom == 0:
                outStr += " O  " + " WAT"
            elif atom == 1:
                outStr += " H1 " + " WAT"
            elif atom == 2:
                outStr += " H2 " + " WAT"
            outStr += str(molecule)
            while len(outStr) < 26:
                outStr = outStr[:20] + " " + outStr[20:]
            outStr += "    "
            outStr += "{:8.3f}".format(data[molecule][atom][0])
            outStr += "{:8.3f}".format(data[molecule][atom][1])
            outStr += "{:8.3f}".format(data[molecule][atom][2])
            outStr += "  1.00" + "  0.00"
            outStr += "           "
            if atom == 0:
                outStr += " O  "
            elif atom == 1:
                outStr += " H  "
            elif atom == 2:
                outStr += " H  "
            outStr += "\n"
            inFile.write(outStr)


# Writes data to PDB file style = HETATOM
def writeDataPDBHETATM(data, inFile):
    iterator = 0
    for molecule in range(len(data)):
        for atom in range(3):
            iterator += 1
            outStr = "HETATM"
            outStr += str(iterator)
            while len(outStr) < 11:
                outStr = outStr[:6] + " " + outStr[6:]
            outStr += " "
            if atom == 0:
                outStr += " O  " + " WAT"
            elif atom == 1:
                outStr += " H1 " + " WAT"
            elif atom == 2:
                outStr += " H2 " + " WAT"
            outStr += str(molecule)
            while len(outStr) < 26:
                outStr = outStr[:20] + " " + outStr[20:]
            outStr += "    "
            outStr += "{:8.3f}".format(data[molecule][atom][0])
            outStr += "{:8.3f}".format(data[molecule][atom][1])
            outStr += "{:8.3f}".format(data[molecule][atom][2])
            outStr += "  1.00" + "  0.00"
            outStr += "           "
            if atom == 0:
                outStr += " O  "
            elif atom == 1:
                outStr += " H  "
            elif atom == 2:
                outStr += " H  "
            outStr += "\n"
            inFile.write(outStr)


# runs program
def testRun(inFile, err, outFile):
    print "Running Test Version of Program..."


# this is the parent runner for the program
def runPgm(inFile, err):
    print "Running Program..."
    data, pdbType = readFile(inFile)
    newData = [ [ [ 0 for i in range(3) ] for j in range(3) ] for k in range(len(data)) ]
    print "Reordering Molecules..."
    for i in range(len(data)):
        positions = getTetrahedronPositions(data[i])
        zeroedMol = newRandOrientation(positions)
        newMol = resetOrientation( data[i][0], zeroedMol )
        newData[i] = newMol
    print "Molecules Reordered"
    connectedMolecules = getNeighbors(newData)  # -1 index = not neighboring
    finalData = newData
    for i in range(len(connectedMolecules)):
        # print "check defects"
        isFine = isDefectiveCheck(err, connectedMolecules, finalData, i)
        # print "isFINE CHECK", isFine
        if not isFine:
            # print "fixing defects"
            while not isFine:
                # print "RerunMol"
                isFine, finalData = rerunMolAndNeighbors(err, connectedMolecules, finalData, i)
                # print "rerunDone"
    writeDataPDB(finalData, pdbType)
    # printData(newData)


badArgs = checkArgs(pdbIN, maxErr, pdbOUT) # stop in case of bad argument

# check input args
if not badArgs: # stop in case of bad argument
    print "Good Arguments, Initializing Reorientiation with", maxErr, "maximum defects"
    # testRun(pdbIN, maxErr, pdbOUT)
    runPgm(pdbIN, maxErr)
elif badArgs:
    print "Bad Arguments, Quitting..."