#!/usr/bin/python

### Author: Gentry Smith, Oklahoma State University
### Created: July 31, 2017, 12PM
### Last Edited: July 31, 2017

### takes file arg with format [ [energy] [pdb_name] ], alters to [ [energy]  [torsion] ], and creates copy with
###     [ [relative energy]  [torsion] ].

import sys


def IOValidator():
    isValid = False
    try:
        inFile = sys.argv[1]
        isValid = True
    except IOError:
        print("Input arg is not a file.\nQuitting...")
        exit()
    return isValid


def GetFileData():
    inData = []
    inFile = open(sys.argv[1], 'r')
    iter = 0
    for line in inFile:
        inLine = line.split()
        inData.append(float(inLine[0]))
        iter = iter + 1
    inFile.close()
    return inData


def Relativize(energies):
    minimum = min(energies)
    # print("Relativize: minimum="+str(minimum))
    newEnergies = []
    for i in range(len(energies)):
        # print("Relativize: index="+str(i))
        # print("Relativize: energy="+str(energies[i]))
        newMin = (float(energies[i]) - float(minimum))
        # print("Relativize: newMin="+str(newMin))
        newEnergies.append((newMin))
        # print("Relativize: newEnergies="+str(newEnergies))
    return newEnergies


def UnifiedScale(energies):
    # print("unifying scale...")
    maxi = max(energies)
    # print("Unify: max=" + str(maxi))
    newEnergies = []
    for i in range(len(energies)):
        # print("Unify: energy=" + str(energies[i]))
        newEner = (float(energies[i]) / maxi)
        # print("Unify: scaled energy=" + str(newEner))
        newEnergies.append(newEner)
    return newEnergies


def CriticalHit(energies, torsions):
    isIncreasing = True
    crits = []
    tors = []
    prev = 0
    for i in range(len(energies)):
        if ( energies[i] == 0 ):
            crits.append(energies[i])
            tors.append(torsions[i])
        if ( (isIncreasing) & (energies[i] < prev) ) or ( (not isIncreasing) & (energies[i] > prev) ):
            crits.append(energies[i-1])
            tors.append(torsions[i-1])
            isIncreasing = not isIncreasing
        prev = float(energies[i])
    returnThing = [crits, tors]
    return returnThing


def MakeFile(energies, torsions, fileName):
    outFile = open(fileName, 'w')
    for i in range(len(energies)):
        strOut = ('{:.11e}'.format(energies[i]) + " " + str(torsions[i]) + "\n")
        outFile.write(strOut)
    outFile.close()


def Runner():
    if IOValidator():
        energies = GetFileData()
        torsions = [180]
        i = 185
        while i != 180:
            if i == 360:
                i = 0
            torsions.append(i)
            i = i + 5
        MakeFile(energies, torsions, 'abs_energ.txt')
        relativeEnergies = Relativize(energies)
        MakeFile(relativeEnergies, torsions, 'rel_energ.txt')
        MakeFile(UnifiedScale(relativeEnergies), torsions, 'uni_energ.txt')
        crits = CriticalHit(relativeEnergies, torsions)
        MakeFile(crits[0], crits[1], 'crit_pts.txt')


Runner()