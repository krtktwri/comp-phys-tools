from myConstants import *
# from myPotentials import *

import cmath as cm # To handle complex  numbers
import numpy as np # To handle arrays

#Function to take modulus of complex number
def mod(c):
    return np.sqrt(c.real**2 + c.imag**2)


dx, dt = 0.01, 0.005

x = np.arange(0, np.pi, dx)
t = np.arange(0, 50, dt)

L = x[-1] - x[0]                                       # Length domain in consideration

Psi = np.zeros((len(t), len(x)), dtype = 'complex_')   # Unknown vector at each time that I need to solve for 

A = np.zeros((len(t), len(x)), dtype = 'complex_') 
B = np.zeros((len(t), len(x)), dtype = 'complex_') 

#Supplementary Matrices
U = np.zeros((len(t), len(x)), dtype = 'complex_') 
R = np.zeros((len(t), len(x)), dtype = 'complex_')

S = np.ones(len(x)-1, float)  

######## GAUSSIAN POTENTIAL ###############
def V_Gauss(x):

    mu = 1.5
    sigma = 0.01
    return 0.1*(1/((np.pi**0.25)*np.sqrt(sigma)))*np.exp(-(x - mu)**2/(2*sigma**2))

###### FINITE WELL POTENTIAL #############
def V_FiniteWell(x):
    
    a, b = 1/3*L, 2/3*L
   
    if x<a or x>b:
        return 3
    if x>=a and x<=b:
        return 0

######## FINITE STEP POTENTIAL #########
def V_FiniteStep(x):
    
    b = (2/3*L)*0.75
    
    if x>b:
        return 0.05
    if x<=b:
        return 0

######## INFINITE WELL POTENTIAL ##########
def V_InfiniteWell(x):
    a, b = 1/3*L, 2/3*L
    
    if x<a or x>b:
        return 10000000
    if x>=a and x<=b:
        return 0


######## FREE PARTICLE POTENTIAL ########    
def V_Free(x):    
    return 0

######## HARMONIC OSCILLATOR POTENTIAL ###############
def V_Harmonic(x):
    w = 1
    return 0.5*m*x*w**2


# Choose what potential you want to evolve here (or define your own)
def V(x):
    
    return V_FiniteStep(x)


# Initial Wave function
mu, sigma = 1, 0.05
p0 = np.sqrt(0.1*m*0.511)                      #average momentum of an electron in chosen units

for i in range(0, len(x)):
    Psi[0, i] = ((1/(np.pi**(1/4)*np.sqrt(sigma)))*np.exp((-(x[i] - mu)**2)/(2*sigma**2))*cm.exp(1j*p0*x[i]/hbar))/10
    
# Potential
p = []
for i in range(len(x)):
    p.append(V(x[i]))    
    
    
    
    
