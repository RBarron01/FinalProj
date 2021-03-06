#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:35:04 2018

@author: Jodie McLennan and Rachel Barron
"""
# tracks the path of the uncharged particles
import numpy as np
import random

lght = 2.99e8 #m/s

def path0(veloc,incRad): 
    # veloc is total velocity magnitude   
    loren = 1/np.sqrt(1-veloc**2/lght**2)
    time = np.arange(0,.0000000005,10**-10)
    pos0X = random.random()*.000001
    pos0Y = random.random()*.000001
    # particles
    posRelX = pos0X + loren*veloc*time*np.cos(incRad)
    posRelY = pos0Y + loren*veloc*time*np.sin(incRad)
    # anti-particles
    antiposRelX = pos0X - loren*veloc*time*np.cos(incRad)
    antiposRelY = pos0Y - loren*veloc*time*np.sin(incRad)
    pos = [posRelX,posRelY]
    antipos = [antiposRelX,antiposRelY]
    return [pos,antipos]
