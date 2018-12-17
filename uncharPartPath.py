#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:35:04 2018

@author: rachelbarron
"""
import numpy as np
import random

lght = 2.99e8

def path0(veloc,incRad): # tracks the path of an uncharged particle
    # veloc is total veloc
    
    loren = 1/np.sqrt(1-veloc**2/lght**2)
    time = np.arange(0,.0000000005,10**-10)
    pos0X = random.random()*.000001
    pos0Y = random.random()*.000001
    posRelX = pos0X + loren*veloc*time*np.cos(incRad)
    posRelY = pos0Y + loren*veloc*time*np.sin(incRad)
    antiposRelX = pos0X - loren*veloc*time*np.cos(incRad)
    antiposRelY = pos0Y - loren*veloc*time*np.sin(incRad)
    pos = [posRelX,posRelY]
    antipos = [antiposRelX,antiposRelY]
    return [pos,antipos]