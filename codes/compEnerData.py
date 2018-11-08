#!/usr/bin/python

### Author: Gentry Smith, Oklahoma State University
### Created: July 31, 2017, 3PM
### Last Edited: August 1, 2017

### Takes data created by teatAbsEnergies and compares values via additive and multiplicative comparison
###     with abs or rel data. Math in terms of File 2 sub/div File 1.

# sys.argv[1] = file 1, working directory here.
# sys.argv[2] = file 2, compared with file 1.


import sys
import numpy
import math

def IOValidator():
    isValid1 = False
    isValid2 = False
    try:
        inFile1 = open(sys.argv[1])
        isValid1 = True
    except IOError:
        print("Arg File 1 is invalid.")
        isValid1 = False
    try:
        inFile1 = open(sys.argv[2])
        isValid2 = True
    except IOError:
        print("Arg File 2 is invalid.")
        isValid2 = False
    if (isValid1 & isValid2 & (sys.argv[1] != sys.argv[2]) ):
        print('Valid Args. Running...')
        return True
    else:
        if (sys.argv[1] == sys.argv[2]):
            print ('args are indentical. Skipping...')
        else:
            print("Invalid args. Quitting...")
        exit()


def ExtractData(data):
    inFile = open(data, 'r')
    inData = []
    inTorsions = []
    # print('Extracting Data...')
    for line in inFile:
        # print('line=' + str(line))
        # print('line.split()=' + str(line.split()))
        # print('line.split()[1]=' + str(line.split()[1]))
        inData.append(float(line.split()[0]))
        inTorsions.append(int(line.split()[1]))
    # print(str(inTorsions))
    # print('Done.')
    return [inData, inTorsions]


def Comparator(data1, data2, func):
    # func: 0=add, 1=mult
    newData = []
    if func == 0:
        for i in range(len(data2)):
            newData.append(float(data2[i] - data1[i]))
    elif func == 1:
        for i in range(len(data2)):
            try:
                newData.append(float(data2[i] / data1[i]))
            except ZeroDivisionError:
                newData.append(0.0)
    return newData


def WriteFile(data1, data2, tors, compData, comp, sigs):
    # writes data of comparison. Format:
    #   File1 = {file1}
    #   File2 = {file2}
    #   Source: {absolute, relative}
    #   Comparison: {additive, multiplicative}
    #   comp: {min/max/avg/stdev of all comp values}
    #   Raw Data: {includes header of File1, File2, Torsions, Comp defining each column}
    # print("Writing file...")
    # print('File2=' + str((sys.argv[2]).split("/")))
    source = ""
    if str(sys.argv[1])[:3] == "abs":
        source = "absolute"
    elif str(sys.argv[1])[:3] == "rel":
        source = "relative"
    elif str(sys.argv[1])[:3] == "uni":
        source = "unified relative scale"
    else:
        print(str(sys.argv[1])[:2])
    comparison = ""
    if comp == 0:
        comparison = "additive"
    elif comp == 1:
        comparison = "multiplicative"
    headerLines = [0]*10
    headerLines[0] = ('File1 = ' + sys.argv[1] + '\n')
    headerLines[1] = ('File2 = ' + sys.argv[2] + '\n')
    headerLines[2] = ('Source: ' + source + '\n')
    headerLines[3] = ('Comparison: ' + comparison + '\n')
    headerLines[4] = ('Comparison min: ' + str(sigs[0]) + '\n')
    headerLines[5] = ('Comparison max: ' + str(sigs[1]) + '\n')
    headerLines[6] = ('Comparison avg: ' + str(sigs[2]) + '\n')
    headerLines[7] = ('Comparison stdev: ' + str(sigs[3]) + '\n')
    headerLines[8] = ('Raw Data:' + '\n')
    f1ColSize = len(str(data1[0]))
    f2ColSize = len(str(data2[0]))
    headerLines[9] = ('File1'.ljust(18) + 'File2'.ljust(18) + 'Tors'.ljust(5) + 'Comp'.ljust(18) + '\n')
    fileName = (str((sys.argv[2]).split("/")[-2]) + "_" + str(sys.argv[1])[:3] + "_" + comparison + '.txt')
    outFile = open(fileName, 'w')
    for i in range(len(headerLines)):
        outFile.write(str(headerLines[i]))
    for i in range(len(data1)):
        # print('str(tors[i]).ljust(5)=' + str(tors[i]).ljust(5))
        string = (str(data1[i])[:17].ljust(18) + ' ' + str(data2[i])[:17].ljust(18) + str(tors[i]).ljust(5) + str(compData[i])[:17].ljust(18) + '\n')
        outFile.write(string)


def GetCompSigs(data):
    sigs = []
    sigs.append(min(data))
    sigs.append(max(data))
    sigs.append((float(sum(data))/float(len(data))))
    sigs.append(numpy.std(data, axis=0))
    return sigs


def Runner():
    if IOValidator():
        [data1, torsions1] = ExtractData(sys.argv[1])
        [data2, torsions2] = ExtractData(sys.argv[2])
        if (len(data1) == len(data2)) & (len(torsions1) == len(torsions2)):
            aData = Comparator(data1, data2, 0)
            aSigs = GetCompSigs(aData)
            WriteFile(data1,data2, torsions1, aData, 0, aSigs)
            mData = Comparator(data1, data2, 1)
            mSigs = GetCompSigs(mData)
            WriteFile(data1, data2, torsions1, mData, 1, mSigs)
            print('Complete.')
Runner()
