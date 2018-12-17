#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:23:04 2018

@author: rachelbarron
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from velocFunc import velocSol
import random

def path(energy,m,q):
    velocX = velocSol(3.5,energy)*np.cos(.01) # m/s
    velocY = velocSol(3.5,energy)*np.sin(.01)
    velocArr = np.array([velocX,velocY,0]) # m/s
    posX = random.random()*.000001
    posY = random.random()*.000001
    pos = np.array([posX,posY,0]) # m
    time = np.arange(0,.000000002,10**-10) # s
    magField = np.array([0,0,3.8]) # T
    m = m*1.79e-30
    veloc_pos = np.array([velocArr[0],velocArr[1],velocArr[2],pos[0],pos[1],pos[2]]) # initial conditions
    veloc_antipos = np.array([-velocArr[0],-velocArr[1],-velocArr[2],pos[0],pos[1],pos[2]]) # initial conditions
    
    def dpvdt (veloc_pos,time,m,q,magfield): # veloc_pos = vx,vy,vz,x,y,z
        # calculate cross product
        veloc = np.array([veloc_pos[0],veloc_pos[1],veloc_pos[2]])
        veloc_force = np.cross(veloc,magfield)
        ####system of differential equations####
        # acceleration m^2/s
        dvxdt = (q/m)*(veloc_force[0])
        dvydt = (q/m)*(veloc_force[1])
        dvzdt = (q/m)*(veloc_force[2])
        # velocity m/s
        dxdt = veloc_pos[0]
        dydt = veloc_pos[1]
        dzdt = veloc_pos[2]
        
        return [dvxdt,dvydt,dvzdt,dxdt,dydt,dzdt]
    
    solution = odeint(dpvdt,veloc_pos,time,args = (m,q,magField),mxstep = 50000000)
    solution2 = odeint(dpvdt,veloc_antipos,time,args = (m,-q,magField),mxstep = 50000000)
    return [solution,solution2]
