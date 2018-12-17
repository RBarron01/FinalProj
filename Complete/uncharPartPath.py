#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:35:04 2018

@author: rachelbarron
"""
import numpy as np
import random

lght = 2.99e8 # m/s

def path0(veloc,incRad): # tracks the path of an uncharged particle
    # This is only for uncharged particles. Takes scalar inputs
    
    loren = 1/np.sqrt(1-veloc**2/lght**2) # lorentz factor
    time = np.arange(0,.0000000005,10**-10) # s
    # next 2 lines: becuase the beam cannot be modeled as a point
    pos0X = random.random()*.000001 
    pos0Y = random.random()*.000001
    # list of positions over time step
    posRelX = pos0X + loren*veloc*time*np.cos(incRad)
    posRelY = pos0Y + loren*veloc*time*np.sin(incRad)
        # and the anti-particles:
    antiposRelX = pos0X - loren*veloc*time*np.cos(incRad)
    antiposRelY = pos0Y - loren*veloc*time*np.sin(incRad)
    pos = [posRelX,posRelY]
    antipos = [antiposRelX,antiposRelY]
    return [pos,antipos]
