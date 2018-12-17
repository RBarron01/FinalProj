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

partList.append(['electron',.511])
partList.append(['muonNeg',105.6])
partList.append(['tauNeg', 1777])
partList.append(['proton',938.3])
partList.append(['neutron',939.6])
partList.append(['delta',1232]) # have 4 deltas so totally random choice between those 4
partList.append(['lambda0',1115.7]) # careful, lambda is a keyword
partList.append(['sigmaPos', 1189.4])
partList.append(['sigmaNeut',1192.5])
partList.append(['sigmaNeg',1197.4])
partList.append(['xiNeut',1315]) 
partList.append(['xiNeg',1321])
partList.append(['omegaNeg',1672])
# mesons
partList.append(['pionChar',139.6])# has antiparticle so will produce +- at the same time
partList.append(['pionNeut',135])
partList.append(['kaonChar',493.7]) # + with - antiparticle
partList.append(['kaonNeut',497.7]) 
partList.append(['eta',547.8])
partList.append(['etaPrime',957.6])

partArray = np.array([0.0]*2 for i in range(len(partList)))
for i in range(len(partList)):
    partSub = partList[i]
    names.append(partSub[0])
    energies.append(partSub[1])
someasdf = [names,energies]

# Create a .txt file with D0-Mesons
partData = open('particleData.txt','w')
partData.write("Name    Energy (eV) \n")
[partData.write("%s    %f\n" % (names[i],energies[i])) for i in range(len(names))]
partData.close()