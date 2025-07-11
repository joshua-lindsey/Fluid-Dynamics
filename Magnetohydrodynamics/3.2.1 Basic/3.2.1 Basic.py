# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 20:09:23 2018

@author: joshu
"""


import matplotlib.pyplot as plt
import numpy as np

#Global Variables
g = 9.81
H = 1000
x = 1
ho = 1
V0 = 1
Va = 3
Bz = 1000


# Functions
def omega_R(k):
    omega = (g*H*k**2 + (Va/H)**2 ) ** (1/2)
    return omega
            
def height(k,t):
    coef = (k * H) / omega_R(k)
    alpha = (k * x - omega_R(k) * t) # (np.pi / 180)
    return coef * np.cos(alpha)

def Bfield(k,t):
    coef = (-Bz * V0) / (H * omega_R(k) )
    alpha = (k * x - omega_R(k) * t )
    return coef * np.sin(alpha)


#################################################
#k = np.arange(0.0, 100.0, 1) 
k=.001
t = np.arange(0.0, 100.0, 1) #min, max, step

##################################################
#Plot for velocity

fig, ax = plt.subplots()
ax.plot(t, Bfield(k,t))
ax.set(xlabel='time (s)', ylabel='B_field (T)', title='B_field versus time')
ax.grid()
plt.show()


################################################
#Plot for height

fig, ax = plt.subplots()
ax.plot(t, height(k,t))
ax.set(xlabel='time (s)', ylabel='height (m)',title='Height versus time')
ax.grid()
plt.show()

################################################
#Plot for omega_R

fig, ax = plt.subplots()
ax.plot(k, omega_R(k))
ax.set(xlabel='k (rad/m)', ylabel='omega_R (rad/s)', title='Dispersion Relation')
ax.grid()
plt.show()
