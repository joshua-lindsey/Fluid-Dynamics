'''
author: Josh Lindsey
date: Oct. 216, 2018
'''

import math
import numpy as np


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

G = 9.81 #m/s^2
H = 1000

###############################################################################
# Functions to calculate coefficients
#kx is wavenumber is x direction

def C1():
    return(8)
    
def C2(eta, kx):
    c2 = 8*eta*kx**2
    return(c2)
    
def C3(eta, kx, Va):
    term1 = G*H*kx**2
    term2 = (Va/H)**2
    term3 = (eta*kx**2)**2
    c3 = 2*(term1 + term2 + term3)
    return(c3)
  
def C4(eta, kx, Va):
    c4 = eta*(kx*Va/H)**2
    return(c4)

###############################################################################
#Function to solve cubic formula, solve function, with corresponding coefficients     
def cubic_solver(eta, kx, Va):
    return solve(C1(), C2(eta, kx), C3(eta, kx, Va), C4(eta, kx, Va))

def real_array_sol(nump_array):
    Reals =[]
    np_array = nump_array
    for i in range(0,3):
        real_elem = np.real(np_array[i])
        complex_elem = np.imag(np_array[i])
        
        if complex_elem == 0.0 :
            Reals.append(real_elem)
        else:
            Reals.append(0)
        return Reals
    
###############################################################################
# empty dicts, cubic returns 3 solutions  
        
sol1 ={}
sol2 ={}
sol3= {}
 
for i in range(0,10):
    for j in range(0,3):
        for k in range(1000, 1005):
            key = i,j,k
            array_of_sol = cubic_solver(i,j,k)
            R_array_of_sol = real_array_sol(array_of_sol)
            sol1.update({key : R_array_of_sol[0] })
    
pos_sol1 = {}

for key,value in sol1.items():
    if value > 0:
        pos_sol1[key] = value
    else:
        pos_sol1[key] = 'no'
        
print(pos_sol1)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   