#!/usr/bin/python

### Author: Gentry Smith, Oklahoma State University
### Created: March 11, 2018 4AM
### Last Edited: March 11, 2018

### Edits a prepared oopse rose water simulation and creates the desired surface by
### generating a new water.md file labeled newWater.md.
### This adjusts distance and quantity for one or two atom types.

# input:
# oopseSurfaceMaker.py [quantity of atoms] [distance between atoms] [quantity of atom types]

import sys

def NoArgRunner():
    print('Easy input:\npython oopseSurfaceMaker.py [quantity of atoms] [distance between atoms] [quantity of atom types] ')
    # print('Welcome to the OOPSE rose surface generator! This adjusts the surface to fit your desired system.')
    print("\nThis script allows up to two types of atoms on the surface.\n")
    # print('Entering non-int values will quit the program.\n')
    # try:
    #     numAtoms = int(raw_input("input the int value of the number of atoms:\n"))
    #     numTypes = int(raw_input("input the int value of the number of atom types:\n"))
    # except ValueError:
    #     print("\n\n\nNot an int. Quitting. . . \n\n\n")
    #     quit()
    # for i in range(1, numTypes + 1):
    #     raw_input('Looking at atom type ' + str(i) + '.\n ')

def ExtractData(data):
    inFile = open(data, 'r')
    dataCheck = 0
    afterData = []
    # print('Extracting Data...')
    for line in inFile:
#	print("'" + line + "'")
#	print(line[:9])
        # print('line=' + str(line))
        # print('line.split()=' + str(line.split()))
        # print('line.split()[1]=' + str(line.split()[1]))
        if (str(line[:9]) == 'molecule{') and (dataCheck <2):
            #print('FOUND: ' + line)
            dataCheck += 1
        elif dataCheck >= 2:
 #           print(line)
	    afterData.append(str(line))

    # print(str(inTorsions))
    # print('Done.')
    return afterData


def GetPosData(data):
    numAtoms = data[0]
    distAtoms = data[1]
    centerIndex = -1
    atomPosData = [0]
    dist = 0.0
    for i in range(1,numAtoms):
        dist = dist + distAtoms
        atomPosData.append(float(dist))
    halfDist = dist / 2
    for i in range(len(atomPosData)):
        atomPosData[i] = (atomPosData[i] - halfDist)
    return atomPosData

    ### creates this data string for each atom:
    # atom[i]{
    #   type = "X";
    #   position( P, 0.0, 0.0 );
    # }
def GetPosDataString(PosData, numTypes):

    # init vars
    posString = []
    typeIter = 0

    # make type strings
    t1 = "HEAVY0"
    t2 = "HEAVY1"
    # t3 = "HEAVY2"
    atomTypes = [ t1, t2 ]

    for i in range(len(PosData)):
        atomString = '  atom[' + str(i) + ']{\n    type = "' \
                     + str(atomTypes[typeIter]) + '";\n    position( ' \
                     + str(PosData[i]) + ', 0.0, 0.0 );\n  }'
        if typeIter == (numTypes - 1):
            typeIter = 0
        else:
            typeIter = typeIter + 1
        posString.append(atomString)

    # make trailing data
    iterStr = '0'
    for i in range(1, len(PosData)):
        iterStr = (str(iterStr) + ', ' + str(i))
    posString.append("\n   rigidBody[0]{\n    members(" + str(iterStr) + ");\n  }\n")
    posString.append("}\n\nmolecule{\n")


    return posString



def WriteWaterMDData(PosString):
    # define standard string data
    headString = '#ifndef _WATER_MD_\n#define _WATER_MD_\n\nmolecule{\n  name = "PART_WALL";\n\n'
    tailString = ExtractData('water.md')
    newWater = open('newWater.md', 'w')
    finalWrite = [headString]
    for each in PosString:
        finalWrite.append(each + '\n')
    for line in tailString:
        finalWrite.append(line)
    for i in finalWrite:
        #if i[-2:] != '\n':
         #   i = i + '\n'
        newWater.write(i)
    newWater.close()


def ArgRunner():
    # ingest input data
    try:
        numAtoms = int(sys.argv[1])
        distAtoms = float(sys.argv[2])
        numTypes = int(sys.argv[3])
        # atomTypeData = []

        # for i in range(1, numTypes + 1):
        #     newData = [ int( sys.argv[ (3+(i*1)) ] ), int( sys.argv[ (3+(i*2)) ] ) ]
        #     atomTypeData.append(newData)
    except IndexError:
        print('\n\nImproper input format. Quitting...\n\n')
        exit()

    # Figure out surface atom positions and write new water.md file
    atomPositions = GetPosData([ numAtoms, distAtoms, numTypes ])
    posString = GetPosDataString(atomPositions, numTypes)
    WriteWaterMDData(posString)

    for i in range(numTypes):
        print('To change the HEAVY' + str(i) + ' charge, enter:\n')
        print("sed 's\// HEAVY0" + str(i) + "          0.0\HEAVY0         X.X\g' DUFF2.frc >> newDuff.frc; mv newDuff.frc DUFF2.frc\n\n")


def Runner():
    if len(sys.argv) == 1:
        NoArgRunner()
    if len(sys.argv) > 1:
        ArgRunner()

Runner()
