#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:12:16 2018

@author: rachelbarron
"""

import numpy as np
from scipy.optimize import bisect

def velocSol(entEnBeam,resMEn):
    # determines the velocity of a particle
    lght = 2.99e8 # m/s
    startGuess = 1000 # m/s
    endGuess = lght/10 # m/s
    resMEn = np.abs(resMEn) # makes sure mass is positive
    def func(veloc):    
        return np.sqrt(entEnBeam/(((1/np.sqrt(1-veloc**2/lght**2))-1)*resMEn)) - veloc
    return bisect(func,startGuess,endGuess,maxiter = 10000000) 
    # increase iterations for wide velocity range 

