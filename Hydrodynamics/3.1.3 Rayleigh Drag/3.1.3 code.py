# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 22:33:48 2018

@author: joshu
"""
'''
Must plot height and omega real seperate from omega imaginary under assumption 
that omega drag is a constant.
When plotting omega_I one can see that the plot is linear and decreasing.
'''

import matplotlib.pyplot as plt
import numpy as np

#Global Variables
g = 9.81
H = 1000
x = 1
ho = 0.1
V_0 = 0.5
omega_drag = 0.1

# Functions
def omega_R(k):
    omega = np.sqrt(g*H*k**2 - (omega_drag /2)**2) 
    return omega

def omega_I():
    omega = - omega_drag/ 2
    return omega
            
def height(k, t):
    coef = (V_0 / g * k) * np.exp(omega_I() * t)
    alpha = (k * x - omega_R(k) * t) #* (np.pi / 180)
    return coef * (omega_R(k) * np.cos(alpha) + omega_I() * np.sin(alpha))



#################################################
k = np.arange(0.0, 100.0, 1)     
t = np.arange(0.0, 100.0, 1) #min, max, step

################################################
#Plot for height

fig, ax = plt.subplots()
ax.plot(t, height(0.25/2*np.pi , t))
ax.set(xlabel='time (s)', ylabel='height (m)',title='Height versus time')
ax.grid()
plt.show()

################################################
#Plot for omega_R

fig, ax = plt.subplots()
ax.plot(k, omega_R(k))
ax.set(xlabel='k', ylabel='omega_R', title='Dispersion Relation')
ax.grid()
plt.show()

################################################
#Plot for omega_I
'''
omega_drag = np.arange(0.0, 100.0, 1)

fig, ax = plt.subplots()
ax.plot(k, omega_I())
ax.set(xlabel='k', ylabel='omega_I', title='Dispersion Relation')
ax.grid()
plt.show()
'''