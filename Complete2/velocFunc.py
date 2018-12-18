#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:12:16 2018

@author: Jodie McLennan and Rachel Barron
"""
# determines the velocity of a particle
import numpy as np
from scipy.optimize import bisect
def velocSol(entEnBeam,resMEn):    
    lght = 2.99e8 # m/s
    startGuess = 1000 # m/s
    endGuess = lght/10 # m/s
    resMEn = np.abs(resMEn) # makes sure mass is positive
    def func(veloc):    
        return np.sqrt(entEnBeam/(((1/np.sqrt(1-veloc**2/lght**2))-1)*resMEn)) - veloc
    # increase iterations in bisect for wide velocity range 
    return bisect(func,startGuess,endGuess,maxiter = 10000000) 
