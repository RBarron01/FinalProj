# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 13:47:33 2018

@author: jmcle
"""
import numpy as np
me = .511 #MeV
mmu = 105.7
mp = 139.6
mk = 493.7
mpr = 938.3

print('input energy of beam(1-1e6)')
Ee = float(input())

positione = np.zeros(100)
positionp = np.zeros(100)
positionp[:] = 50
Erange = np.arange(1,1e6,1)
anglerange = np.arange(80,2)
velocitye = 10e6
Ethreshu = 2*(mmu)**2/me-me

for i in np.arange(0,len(positione)):
    if Ee < Ethreshu:
        positione[i] = positione[i-1] + velocitye
        if positione[i] == positionp[i]:
            man ='man'

            
            
            

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import bisect



# constants from article c71092
bendRad = 100 # m
beamPow = 5.3 #pv megaWatts
circum = 400*(np.pi + 1) # circumference of collider
deltV0 = .025 # experimental and theoretical constant
entEn = 15000 # MeV kinetic energy of e+- when it enters the collider

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
veloc = np.arange(2.22e5,2.3e5,1) # range narrowed through experimentation
def func(veloc):    
    return np.sqrt(entEnBeam/((1/np.sqrt(1-veloc**2/lght**2))-1)/elResEn) - veloc

funcOut = np.array([func(veloc[i]) for i in range(len(veloc))])
plt.plot(veloc,funcOut)
solVeloc = bisect(func,226335,226336,maxiter = 10000)



# calculate momentum
loren = 1/np.sqrt(1-solVeloc**2/lght**2)
momPartMag = loren*elResEn*solVeloc
# the two particles have opp x momentums but similar y 
momElX = momPartMag*np.cos(incRad)
momElY = momPartMag*np.sin(incRad)

momPosX = momPartMag*np.cos(-incRad)
momPosY = momPartMag*np.sin(incRad)


'''
Next steps:
    calculate probabilities from threshold energies
        to follow up, does a particle take more from the e+ or e-?
    calculate path of paticles from magnetic field 
        Includes charge on respective particles





buNum = 6 # number of bunches @Fermilab tevatron number of bunches in 1 dir
bunLen = .06 # m (can range from a few to 50 cm)
sigX = .1 # cm gaussian rms radii wrt x: for now assume all radii to be equal
sigY = .01 # cm gaussian rms radii wrt y


lght = 2.99e10 # m/s

#ken = .5*m*omega**2
enPart = 2#/6.242e9 # GeV

#potEn = 88e3*enPart**4/(lght**2*bendRad)
enLoss = (88/lght**2)
freq = 10**6 # sec^-1 rev freq
iBun = 10**-3 # mA/bunch at 2GeV

massPart = 9.11e-31
gammaPart = enPart/massPart # we think gamma = c^2 = E/mass_e-

radPart = 2.82e-13 # m radius of particles (electron and positron)
 
betaX = np.linspace(.1,1,100) # amplitude functions at interaction point

crossA = 4*np.pi*sigX*sigY
# Think that these vals are wrong... check interp of nums
nElMax = gammaPart*2*np.pi*sigX*(sigX + sigY) * .06/(radPart*betaX)
nPosMax = gammaPart*2*np.pi*sigX*(sigX + sigY) * .06/(radPart*betaX)

eChar = 1.6027e-19 # Coulombs
k = 1/(4*np.pi*eChar**2)

lum = .5e32 # cm^-2*s^-1#nElMax*nPosMax*buNum*freq/(4*np.pi*sigX*sigY)
evRat = lum*crossA # event rate N

iDens = nPosMax*freq*buNum'''
