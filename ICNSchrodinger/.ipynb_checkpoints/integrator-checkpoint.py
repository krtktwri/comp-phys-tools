import cmath as cm
import numpy as np 

from myConstants import *
from initializing import *

def ICNevolve():
    for n in range(0, len(t)-2):              #Index n iterates of the time variables
        Psi[n, 0], Psi[n, -1] = 0, 0          #Fixed Boundary conditions (can be interpreted as infinite potential walls on each sides)

        for i in range(1, len(x)-2):          #Index i iterates over the space variables 
            A[n, i] = -2  + (4j*m*dx**2)/(hbar*dt) - (2*m*dx**2)/(hbar**2)*V(x[i])
            B[n, i] = -Psi[n, i+1] - Psi[n, i-1] + Psi[n, i] * (2 + ((4j*m*dx**2)/(hbar*dt)) + ((2*m*dx**2)/(hbar**2))*V(x[i])) 

        #Supplementary Matrices required for 
        U[n, 1] = 1/A[n, 1]
        R[n, 1] = B[n, 1] * U[n, 1]

        #Foward Sweep
        for i in range(1, len(x)-2):
            U[n, i] = 1/(A[n, i] - U[n, i-1])
            R[n, i] = (B[n, i] - R[n, i-1])*U[n, i]

        N = len(x)-1
        i = N-1

        Psi[n+1, N] = R[n, N] 

        #Backward Sweep
        while i>=1:
            Psi[n+1, i] = R[n, i] - U[n, i]*Psi[n+1, i+1]
            i -= 1
    
    return Psi

############################################################

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os


def animate():

    # Potential
    po = []
    for i in range(len(x)):
        po.append(V(x[i]))    

    fig = plt.figure()
    plt.xlabel('Position [x]')
    plt.ylabel('Probablity Distribution [U]')
    plt.title("Free Particle Evolution")
    plt.grid()

    plts = []             
    # plt.hold()
    for i in range(0, int(len(t)), 5):
        p, = plt.plot(x, mod(Psi[i,:]**2), 'r')   
        q, = plt.plot(x, po, 'kx--')      
        plts.append( [p, q] )           

    ani = animation.ArtistAnimation(fig, plts, interval=1, repeat_delay=3000)   # run the animation
    ani.save('anim/animation.mp4', fps=24)
    os.system("ffmpeg -i QMAnimation.mp4 QMAnimation.gif")