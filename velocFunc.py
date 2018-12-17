#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:12:16 2018

@author: rachelbarron
"""

import numpy as np
from scipy.optimize import bisect

def velocSol(entEnBeam,resMEn):
    lght = 2.99e8
    startGuess = 1000
    endGuess = lght/10
    resMEn = np.abs(resMEn)
    def func(veloc):    
        return np.sqrt(entEnBeam/(((1/np.sqrt(1-veloc**2/lght**2))-1)*resMEn)) - veloc
    return bisect(func,startGuess,endGuess,maxiter = 10000000) 

