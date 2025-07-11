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

#Do not touch code below
###############################################################################
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


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0
# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0
# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)
###############################################################################
###############################################################################

#my code is going to start here
#Function of C represent coefficients of the cubic equation
#my code is going to start here
#Function of C represent coefficients of the cubic equation

G = 24.5 #m/s^2 Jupiter
H = 5000

###############################################################################
# Functions to calculate coefficients

def C1():
    c1 = -8
    return(c1)
    
def C2(eta, Fo, k, wv):
    c2 = -8*(wv+8*eta*k**2-8*Fo)
    return(c2)
    
def C3(eta, Fo, k, Va, wv):
    t1 = G*H*k**2
    t2 = (Va/H)**2
    t3 = 3*wv*eta*k**2
    t4 = -3*Fo*wv
    t5 = -3*Fo*eta*k**2
    t6 = wv**2
    t7 = Fo**2
    t8 = eta**2 * k**4
    c3 = (-2)*(t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8)
    return(c3)
  
def C4(eta, Fo, k, Va, wv):
    t1 = (wv + eta *k**2)*(Fo**2 + (Va/H)**2)
    t2 = (wv - Fo) * (G*H*k**2 + eta**2 * k**4)
    t3 = wv**2 * (eta*k**2 -Fo)
    t4 = -2*wv*Fo*eta*k**2
    t5 = 2*G*H*eta*k**4
    t6 = -2*Fo*(wv*eta*k**4 + (Va/H)**2)
    c4 = t1 + t2 + t3 + t4 + t5 + t6
    return(c4)

###############################################################################
#Function to solve cubic formula, solve function, with corresponding coefficients     
def cubic_solver(eta, Fo, k, Va, wv):
    return solve(C1(), C2(eta, Fo, k, wv), C3(eta, Fo, k, Va, wv), C4(eta, Fo, k, Va, wv))

def real_array_sol(nump_array):
    Reals =[]
    np_array = nump_array
    for i in range(0,3):
        real_elem = np.real(np_array[i])
        complex_elem = np.imag(np_array[i])
        #print(real_elem,'real elem')
        #print(complex_elem,'complex elem')
        if (real_elem > 0.0 and complex_elem == 0):
            Reals.append(real_elem)
        else:
            Reals.append(0)
            
    return Reals

# This is the end of the functions    
###############################################################################
###############################################################################
###############################################################################        
###############################################################################
    
#print(cubic_solver(1,1,1,1,1))

# empty dict for three solutions        
sol1 ={}
sol2 ={}
sol3= {}
 

# v1 = eta
# v2 = Fo
# v3 = k
# v4 = Va
# v5 = wv

# Timer start
start1 = timeit.default_timer()


step1 = 50
startpoint1 = 450
endpoint1 = 500

for v1 in arange(0, 100, 25):
    for v2 in arange(0, 250, 50):
        for v3 in arange(0, 200, 25):
            for v4 in arange(0, 100, 20):
                for v5 in arange(100, 250, 50):
                    key = v1, v2, v3, v4, v5
                    array_of_sol = cubic_solver(v1, v2, v3, v4, v5)
                    #print(array_of_sol,'array of sol')
                    R_array_of_sol = real_array_sol(array_of_sol)
                    #print(R_array_of_sol,' real array of sol')
                    #if R_array_of_sol[0] > 0:
                    sol1.update({key : R_array_of_sol[0] })
#print(sol1, ' sol1')
    
pos_sol1 = {}

#The if statement checks for positive values
for key,value in sol1.items():
    if value > 0:
        pos_sol1[key] = value
    else:
        pass

if not pos_sol1:
    print('no solutions')
else:
    print(pos_sol1) 

stop1 = timeit.default_timer()  
print(round( stop1-start1 ,8), 's')

print('boop booop computer done...')    
'''
###############################################################################    
###############################################################################
# block of code unpacks dict into 4 list of x,y,z,vals

cord = list(pos_sol1.keys())
vals = list(pos_sol1.values())

V1 = []
V2 = []
V3 = []
V4 = []
V5 = []

#This unpacks keys into x,y,z
for i in range(len(cord)):
    V1.append(cord[i][0])
    V2.append(cord[i][1])
    V3.append(cord[i][2])
    V4.append(cord[i][3])
    V5.append(cord[i][4])
    
#print(x,y,z,vals)
###############################################################################
start1 = timeit.default_timer()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter(V1,V2,V3, c='g', marker='o')

ax.set_xlabel('X (eta)' , fontsize = 10)
ax.set_ylabel('Y (Fo)' , fontsize = 10)
ax.set_zlabel('Z (k)', fontsize = 10)
ax.dist =12
fig = plt.figure(figsize=(4.,35.))
ax.set_title('Positive Real Roots Solution #1', loc='center')
plt.show()    
    
stop1 = timeit.default_timer()  

#print(round(stop1-start1 , 8), 's')    


'''