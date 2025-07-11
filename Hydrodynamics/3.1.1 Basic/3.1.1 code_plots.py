# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 10:45:28 2018

@author: joshu
"""

import matplotlib.pyplot as plt
import numpy as np

#Global Variables
g = 9.81
H = 1000
x = 1
ho = 0.1
V_0 = 0.1


# Functions
def omega_R(k):
    omega = np.sqrt(g*H) * k
    return omega
            
def height(k,t):
    coef = (k * H * V_0) / omega_R(k)
    alpha = (k * x - omega_R(k) * t) #* (np.pi / 180)
    return coef * np.cos(alpha)

def velocity(k,t):
    alpha = (k * x - omega_R(k) * t )
    a =  V_0 * np.cos(alpha)
    return a


#################################################
#k = np.arange(0.0, 200.0, 1) 
k=2    
t = np.arange(0.0, 200.0, 1) #min, max, step

##################################################
#Plot for velocity

fig, ax = plt.subplots()
ax.plot(t, velocity(k,t))
ax.set(xlabel='time (s)', ylabel='velocity (m/s)', title='Velocity versus time')
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
