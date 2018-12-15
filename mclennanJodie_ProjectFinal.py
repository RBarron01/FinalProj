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
