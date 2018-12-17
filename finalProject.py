#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 12:13:26 2018
@author: Rachel Barron and Jodie McLennan
"""
#import the modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import bisect
import pandas as pd
import partTxtGen #creates txt file
from astropy.table import Table, Column
import random
from scipy.integrate import odeint

#functions we created
from velocFunc import velocSol
from partPath import path
from uncharPartPath import path0



# constants from article c71092
bendRad = 100 # m
beamPow = 5.3 #pv megaWatts
circum = 400*(np.pi + 1) # circumference of collider
deltV0 = .025 # experimental and theoretical constant
entEn = 15000 # MeV kinetic energy of e+- when it enters the collider
bField = [0,0,3.8] # T


lght = 2.99e8 # m/s
elResEn = .511 #MeV
nPart = 1e11 # particles/bunch

radEl = 2.82e-15 # m radius of particles (electron and positron)
massEl = 9.11e-31 # kg

betaE = .05 # m beta func at interaction point
#effCrossA = 5e-39 # m^2 interaction area of beams
incRad = .01 # rad rad vert angle of incindence between e+- beams
beamWid = .014 # m
beamRad = beamWid/2 # m

# Calculations
lumEE = (3/(16*np.pi))*deltV0*bendRad*beamPow/((radEl**2)*massEl*(lght**2)*(elResEn**3)*betaE)
lee = .5e32
gamma1 = entEn/elResEn
gamma = ((3/(16*np.pi))*deltV0*bendRad*beamPow/(radEl**2*massEl*(lght**2)*(lee)*betaE))**(1/3)
effCrossVol = np.pi*beamWid**3*np.cos(incRad)/(2*np.sin(incRad)) # m^3
lenInt = 2*beamWid*np.cos(incRad)/np.sin(incRad) # 

kEn = entEn - elResEn*lght**2 # MeV

# solve for velocity of the particles
entEnBeam = entEn/2
solVeloc = velocSol(entEnBeam,elResEn)


# calculate momentum
loren = 1/np.sqrt(1-solVeloc**2/lght**2) # relativistic lorentz factor
momPartMag = loren*elResEn*solVeloc
# the two particles have opp x momentums but similar y 
momElX = momPartMag*np.cos(incRad)
momElY = momPartMag*np.sin(incRad)

momPosX = momPartMag*np.cos(-incRad)
momPosY = momPartMag*np.sin(-incRad)

# calculate the probabilities from threshold energies
# Particle threshold energy
partList = []
names = []
energies = []

partList.append(['electron',.511])
partList.append(['muonNeg',105.6])
partList.append(['tauNeg', 1777])
partList.append(['proton',938.3])
partList.append(['neutron',939.6])
partList.append(['delta2Pos',1232]) # have 4 deltas so totally random choice between those 4
partList.append(['deltaPos',1232])
partList.append(['deltaNeut',1232])
partList.append(['deltaNeg',1232])
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
partList.append(['kaonNeut',497.7]) # has anti-neutral pair
partList.append(['eta',547.8])
partList.append(['etaPrime',957.6])

echar = -1.602e-19
charges = [echar,echar,echar,-echar,0,-2*echar,-echar,0,echar,0,-echar,0,echar,0,echar,echar,-echar,0,-echar,0,0,0]

partArray = np.array([0.0]*2 for i in range(len(partList)))
for i in range(len(partList)):
    partSub = partList[i]
    names.append(partSub[0])
    energies.append(partSub[1])
antinames = []
[antinames.append('anti-'+i) for i in names]
    
totEnProbs = sum(energies)
deciProbs = []

[deciProbs.append((totEnProbs-i)/totEnProbs) for i in energies]

deciProbs = [np.round(i) for i in deciProbs]
deciProbsTot = sum(deciProbs)

totParts = 100

postii = 0
munegee=0
taunegg=0
randNums = [np.array(random.randint(0,deciProbsTot)) for i in range(totParts)]
partCount = np.array([0.0]*len(deciProbs))
for i in randNums:
    for j in range(len(deciProbs)):       
        if i>=sum(deciProbs[0:j]) and i<=sum(deciProbs[0:j+1]):
            partCount[j]+=1
            break
        elif i<deciProbs[0]:
            partCount[0]+=1
            
count = 0
time = np.linspace(0,1e8,100)
plt.figure()

for i in range(len(partCount)):
    
    if charges[i] == 0:
        xMom1 = -momElX
        yMom1 = -momElY
        xMom2 = -momPosX
        yMom2 = -momPosY
        velocNeut = velocSol(entEnBeam,energies[i])
        posNeutTot = path0(velocNeut,incRad) 
        antiposNeut = posNeutTot[0]
        posNeut = posNeutTot[1]
        plt.plot(posNeut[0],posNeut[1], label = names[i])
        plt.plot(antiposNeut[0],antiposNeut[1], label = antinames[i])
        
    else:
        
        vPart = -velocSol(entEnBeam,energies[i]) # mass in rest energy units
        solutionTot = path(energies[i],np.abs(energies[i]),charges[i])
        solution = solutionTot[0]
        antisolution = solutionTot[1]
        plt.plot(solution[:,3],solution[:,4], label = names[i])
        plt.plot(antisolution[:,3],antisolution[:,4], label = antinames[i])

        
plt.xlabel('x',color = 'b')
plt.ylabel('y',color = 'b')
plt.title('Position of particles from collisions') 
#plt.legend()       
# bar graph
plt.figure()
y_pos = np.arange(len(names))
plt.bar(y_pos, partCount, align='center', alpha=0.5)
plt.xticks(y_pos, names, rotation = 80)
plt.ylabel('Number of Particles')
plt.title('What Particles Came Out of Collision')
plt.show()
   
