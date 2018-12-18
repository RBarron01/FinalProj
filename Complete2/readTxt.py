#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 17:08:04 2018

@author: rachelbarron
"""

import partTxtGen

def readFile():
    names = []
    energies = []
    charges = []
    data = open('particleData.txt','r')
    data.readline() # to skip first line
    echar = -1.602e-19
    for line in data:
        line = line.split()
        for lne in line:
            if line[0] == lne:
                names.append(str(lne))
            elif line[1] == lne:
                energies.append(float(lne))
            else:
                charges.append(float(lne)*echar)
    return [names,energies,charges]