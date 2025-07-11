# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 19:02:46 2018

@author: joshu
"""
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

#Global Variables
g = 24.5
H = 1000
x = 1
ho = 0.1
V0x = 0.1
Fo = 7500
Nu = 1
wd = 1000

# Functions
def omega_R(k):
    omega = np.sqrt(g*H*k**2 - ( (Fo + Nu * k**2 + wd )/2)**2) 
    return omega

def omega_I(k):
    omega = ( Fo - (Nu * k**2 + wd) ) /2
    return omega
            

def height(k,t):
    coef = V0x / ( g*k ) * np.exp(omega_I(k) * t)
    alpha = (k * x - omega_R(k) * t) #* (np.pi / 180)
    return coef *  (omega_R(k) * np.cos(alpha) / omega_I(k) * np.sin(alpha) )

#################################################
k = np.arange(0.0, 100.0, 1) 

#k=2    
t = np.arange(0.0, 100.0, 1) #min, max, step

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

################################################
#Plot for omega_R

fig, ax = plt.subplots()
ax.plot(k, omega_I(k))
ax.set(xlabel='k (rad/m)', ylabel='omega_I (rad/s)', title='Dispersion Relation')
ax.plot([100,0],[0,0],'k')
ax.grid()

plt.fill_between(k, 0, omega_I(k), where=(omega_I(k)) > 0 , color='lightskyblue')

red_patch = mpatches.Patch(color='lightskyblue', label='Instabilities')
plt.legend(handles=[red_patch])


plt.show()























