# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 22:13:29 2018

@author: joshu
"""

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
Nu = 1
wd = 1

# Functions
def omega_R(k):
    omega = np.sqrt(g*H*k**2 - ( (Nu * k**2 + wd )/2)**2) 
    return omega

def omega_I(k):
    omega = -( Nu * k**2 + wd) /2
    return omega
            
#################################################
k = np.arange(0.0, 200.0, 1)   #min, max, step

################################################
#Plot for omega_R

fig, ax = plt.subplots()
ax.plot(k, omega_R(k))
ax.set(xlabel='k (rad/m)', ylabel='omega_R (rad/s)', title='Dispersion Relation')
ax.grid()
plt.show()

################################################
#Plot for omega_I

fig, ax = plt.subplots()
ax.plot(k, omega_I(k))
ax.set(xlabel='k (rad/m)', ylabel='omega_I (rad/s)', title='Dispersion Relation')
ax.grid()
plt.show()