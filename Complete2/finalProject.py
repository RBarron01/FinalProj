#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 12:13:26 2018
@author: Jodie McLennan and Rachel Barron
"""
import numpy as np
import matplotlib.pyplot as plt
import random

from velocFunc import velocSol
from partPath import path
from uncharPartPath import path0
from readTxt import readFile

### accelorator paramaters ###
entEn = 15000 # MeV kinetic energy of e+- when it enters the collider
entEnBeam = entEn/2 # Kinetic energy per particle
incRad = .01 # (rad) verticle angle of incindence between e+- beams

# calculate the probabilities of production from threshold energies
    
[names,energies,charges] = readFile() # particle data lists

# antiparticle names:
antinames = []
[antinames.append('anti-'+i) for i in names]
 
### plot path of resultant particles ###           
plt.figure() # plot e+- collision results
for i in range(len(energies)):    
    if charges[i] == 0:
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

### particle results for 100 incindents ###
totEnProbs = sum(energies)
deciProbs = []

[deciProbs.append((totEnProbs-i)/totEnProbs) for i in energies]

deciProbs = [np.round(i) for i in deciProbs]
deciProbsTot = sum(deciProbs)

totParts = 100
randNums = [np.array(random.randint(0,deciProbsTot)) for i in range(totParts)]
partCount = np.array([0.0]*len(deciProbs))
for i in randNums:
    for j in range(len(deciProbs)):       
        if i>=sum(deciProbs[0:j]) and i<=sum(deciProbs[0:j+1]):
            partCount[j]+=1
            break
        elif i<deciProbs[0]:
            partCount[0]+=1

# bar graph
plt.figure()
y_pos = np.arange(len(names))
plt.bar(y_pos, partCount, align='center', alpha=0.5)
plt.xticks(y_pos, names, rotation = 80)
plt.ylabel('Number of Particles')
plt.title('What Particles Came Out of Collision')
plt.show()
   
