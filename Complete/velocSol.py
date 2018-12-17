#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:12:16 2018

@author: rachelbarron
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect

def velocSol(entEnBeam,elResEn,resMEn):
    lght = 2.99e8
    veloc = np.arange(2.22e5,2.3e5,1) # range narrowed through experimentation
    def func(veloc):    
        return np.sqrt(entEnBeam/((1/np.sqrt(1-veloc**2/lght**2))-1)/resMEn) - veloc
    
    funcOut = np.array([func(veloc[i]) for i in range(len(veloc))])
    plt.plot(veloc,funcOut)
    plt.xlabel('Velocity candidates (m/s)')
    plt.ylabel('Function to determine velocity')
    return bisect(func,226335,226336,maxiter = 10000) 

