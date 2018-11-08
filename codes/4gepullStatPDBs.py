#!/usr/bin/python

### Author: Gentry Smith, Oklahoma State University
### Created: August 7, 2017, 3PM
### Last Edited: August 7, 2017

### Takes a stationary_points.txt file and will copy .pdb files of the same name from a split_conformers.pdb/ folder
### into a new folder "stationary_conformers"

# This does not use any args and instead relies on the stationary points file being "stationary_points.txt" and the
# conformrs residing in a "split_conformers.pdb/" directory on the same level. It will create the new folder "stationary_conformers"

import os

def IOValidator():
    returnBool = [False,False]
    try:
        file1 = open('stationary_points.txt', 'r')
        file1.close()
        returnBool[0] = True
    except IOError:
        print("Did not find 'stationary_points.txt' file. Quitting...")
        quit()
    try:
        wkdir = os.getcwd()
        file2 = os.chdir('split_conformers.pdb')
        os.chdir(wkdir)
        returnBool[1] = True
    except OSError:
        print("Did not find 'split_conformers.pdb' folder. Quitting...")
        quit()
    if returnBool[0] & returnBool[1]:
        return True
    else:
        return False


def GetPDBs():
    pdbNames = []
    inFile = open('stationary_points.txt', 'r')
    for line in inFile:
        pdbNames.append(line.split()[1])
    return pdbNames


def CopyPDBs(pdbList):
    wkdir = os.getcwd()
    for i in range(len(pdbList)):
        pstring = ( 'cp ' + 'split_conformers.pdb/' + str(pdbList[i]) + ' stationary_conformers/' )
        os.popen(pstring)


def Runner():
    if IOValidator():
        print('Valid Args. Running...')
        pdbList = GetPDBs()
        try:
            os.mkdir('stationary_conformers')
            CopyPDBs(pdbList)
        except OSError:
            print("'stationary_conformers' directory already exists. Erase directory and run again. Quitting...")
            quit()


Runner()

