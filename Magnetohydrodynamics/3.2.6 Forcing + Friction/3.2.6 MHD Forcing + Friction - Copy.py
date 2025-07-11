# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 11:29:21 2018

@author: joshu
"""


import math
import numpy as np
from numpy import arange
import timeit
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#Do not touch Solve Function
###############################################################################
def solve(a, b, c, d):

    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)

        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j
        x13 = np.around(x1, 5, out=None)
        x23 = np.around(x2, 5, out=None)
        x33 = np.around(x3, 5, out=None)
        #print('One Real and two complex roots')
        return np.array([x13, x23, x33])           # Returning One Real Root and two Complex Roots as numpy array.

###############################################################################
# Helper function to return float value of f in Solve Function
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0
# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0
# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)
###############################################################################
#Global Variables
G = 10 #m/s^2 Jupiter
H = 1000

###############################################################################
# Functions to calculate coefficients
#kx is wavenumber is x direction

def C1():
    c1 = 1
    return(c1)
    
def C2(Fo):
    c2 = -Fo
    return(c2)
    
def C3(Fo, kx, Va):
    term1 = Fo**2
    term2 = G*H*kx**2
    term3 = (Va/H)**2
    c3 = (1/4)*(term1 + term2 + term3)
    return(c3)
  
def C4(Fo, kx):
    c4 = -(1/8) * Fo * G * H * kx**2
    return(c4)

###############################################################################
#Function to solve cubic formula, solve function, with corresponding coefficients     
def cubic_solver(Fo, kx, Va):
    return solve(C1(), C2(Fo), C3(Fo, kx, Va), C4(Fo, kx))
#Quick check for dependent variables    
#print(cubic_solver(1,1,1))

#This function finds
def real_array_sol(nump_array):
    Reals =[]
    np_array = nump_array
    for i in range(0,3):
        real_elem = np.real(np_array[i])
        complex_elem = np.imag(np_array[i])
        #print(real_elem, 'real elem', complex_elem, 'complex elem')
        if (real_elem > 0.0 and complex_elem == 0):
            Reals.append(real_elem)
        else:
            Reals.append(0)
    return Reals
  
###############################################################################
#Dump for variables that I'm not sure where else to put

def all_solutions(i_i, i_f, i_stp, j_i, j_f, j_stp, k_i, k_f, k_stp):
    # empty dict for three solutions        
    all_sol ={}

    start = timeit.default_timer()
    
    for i in arange(i_i, i_f, i_stp):
        for j in arange(j_i, j_f, j_stp):
            for k in arange(k_i, k_f, k_stp):
                key = i,j,k
                array_of_sol = cubic_solver(i,j,k)
                #print(array_of_sol,'array of sol')
                R_array_of_sol = real_array_sol(array_of_sol)
                #print(R_array_of_sol,' real array of sol')
                all_sol.update({key : R_array_of_sol[0] })
                
    #print(all_sol, ' sol1')
    stop = timeit.default_timer()  
    print(round(stop-start,8), 's')
    print('boop booop computer done...') 
    return all_sol

def pos_solutions():
    pos_sol1 = {}
    sol = all_solutions()
    #The if statement checks for positive values
    for key,value in sol.items():
        if value > 0:
            pos_sol1[key] = value
        else:
            pass
            
    #print(pos_sol1) 
    return pos_sol1
   

###############################################################################    
###############################################################################
# block of code unpacks dict into 4 list of x,y,z,vals

cord = list(pos_sol1.keys())
vals = list(pos_sol1.values())

x=[]
y=[]
z=[]
v=[]

#This unpacks keys into x,y,z
for i in range(len(cord)):
    x.append(cord[i][0])
    y.append(cord[i][1])
    z.append(cord[i][2])
    
#print(x,y,z,vals)
###############################################################################
start1 = timeit.default_timer()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter(x, y, z, c='b', marker='o')

ax.set_xlabel('X (Fo)' , fontsize = 10)

ax.set_ylabel('Y (k)' , fontsize = 10)
ax.set_zlabel('Z (Va)', fontsize = 10)
ax.dist =12
#ax.zaxis.set_rotate_label(False)
fig = plt.figure(figsize=(4.,35.))
#ax.view_init(45,-45)
ax.set_title('Positive Real Roots Solution #1', loc='center')

plt.show()    
    
stop1 = timeit.default_timer()  

print(round(stop1-start1 , 8), 's')    
