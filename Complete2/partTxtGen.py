#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 08:15:49 2018

@author: rachelbarron
"""
# Just to create a txt file to read in data from
import numpy as np

partList = []
names = []
energies = []
charges = []
echar = -1
## data: particle name, threshold energy, charge (in eV)
    # baryons
partList.append(['electron',.511,echar])
partList.append(['muonNeg',105.6,echar])
partList.append(['tauNeg', 1777,echar])
partList.append(['proton',938.3,-echar])
partList.append(['neutron',939.6,0])
partList.append(['delta2Pos',1232,-2*echar]) # have 4 deltas so totally random choice between those 4
partList.append(['deltaPos',1232,-echar])
partList.append(['deltaNeut',1232,0])
partList.append(['deltaNeg',1232,echar])
partList.append(['lambda0',1115.7,0]) # careful, lambda is a keyword
partList.append(['sigmaPos', 1189.4,-echar])
partList.append(['sigmaNeut',1192.5,0])
partList.append(['sigmaNeg',1197.4,echar])
partList.append(['xiNeut',1315,0]) 
partList.append(['xiNeg',1321,echar])
partList.append(['omegaNeg',1672,echar])
    # mesons
partList.append(['pionChar',139.6,-echar])# has antiparticle so will produce +- at the same time
partList.append(['pionNeut',135,0])
partList.append(['kaonChar',493.7,-echar]) # + with - antiparticle
partList.append(['kaonNeut',497.7,0]) # has anti-neutral pair
partList.append(['eta',547.8,0])
partList.append(['etaPrime',957.6,0])

partArray = np.array([0.0]*2 for i in range(len(partList)))
for i in range(len(partList)):
    partSub = partList[i]
    names.append(partSub[0])
    energies.append(partSub[1])
    charges.append(partSub[2])


# Create a .txt file with D0-Mesons
partData = open('particleData.txt','w')
partData.write("Name    Energy (MeV)    Charge (eV) \n")
[partData.write("%s    %f    %i\n" % (names[i],energies[i],charges[i])) for i in range(len(names))]
partData.close()