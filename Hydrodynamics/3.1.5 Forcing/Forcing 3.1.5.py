# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 18:31:48 2018

@author: joshu
"""

import matplotlib.pyplot as plt
import numpy as np


#Global Variables
g = 9.81
H = 1000
x = 1
ho = 0.1
V0x = 0.1
Fo = 100000


# Functions
def omega_R(k):
    #omega = np.sqrt( g*H*(k**2) )#- (Fo/2)**2 ) 
    #return omega
    return np.sqrt(H* (k**2))

def omega_I(Forcing):
    omega = Forcing / 2 
    return omega
            
            
def height(k,t):
    coef = V0x / ( g*k ) * np.exp(omega_I(k) * t)
    alpha = (k * x - omega_R(k) * t) #* (np.pi / 180)
    return coef *  (omega_R(k) * np.cos(alpha) / omega_I(k) * np.sin(alpha) )

#################################################
k = np.arange(0.0, 25.0, .011) 
#Fo = np.arange(0.0, 2000.0, 1) 
#k=2    
t = np.arange(0.0, 25.0, 1) #min, max, step

################################################
#Plot for height
'''
fig, ax = plt.subplots()
ax.plot(t, height(k,t))
ax.set(xlabel='time (s)', ylabel='height (m)',title='Height versus time')
ax.grid()
plt.show()
'''
################################################
#Plot for omega_R

fig, ax = plt.subplots()
ax.plot(k, omega_R(k))
ax.set(xlabel='k (rad/m)', ylabel='omega_R (rad/s)', title='Dispersion Relation')
ax.grid()
plt.show()

################################################
#Plot for omega_I
'''
fig, ax = plt.subplots()
plt.plot(Fo,omega_I(Fo))
ax.set(xlabel='k (rad/m)', ylabel='omega_I (rad/s)', title='Dispersion Relation')
ax.grid()
plt.show()
'''






