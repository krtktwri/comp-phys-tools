import cmath as cm           # For dealing with complex numbers
import numpy as np           # For dealing with arrays

from initializing import *   # Importing this initialized arrays with the choice of potential distribution and initial wavefunction
from integrator import *     # This module has functions to integrate and animate

import matplotlib.pyplot as plt  # For plotting purposes

#Graph Parameters
plt.rcParams['figure.figsize'] = 15, 10
plt.rc('text', usetex=False)
plt.rcParams.update({'font.size': 20})
plt.rcParams['font.family'] = 'serif'

# Filling the Psi array with solutions from ICN integrator
Psi = ICNevolve()

up, down = (Psi[0].real).max(), (Psi[0].real).min()                   #Adjusting the limit of the graph


# Plotting the initial state of the system
plt.plot(x, (Psi[0]).real, label=r"Wavefunction $\Psi(x)$")
plt.plot(x, mod(Psi[0])**2, label=r"Probablity Distribution $|\Psi(x)|^2$")
plt.plot(x, p, 'k--', label="Potential")
# plt.ylim([down*1.25, up*1.25])
plt.title("Initial State of the System")
plt.xlabel("X")
plt.legend()
plt.grid()
plt.savefig("InitialState.png")


# Plotting time evolution of wavefunction
n = 9              # Number of curves I want to plot
m = int(len(t)/n)  # Interval in terms of time index

#Plotting several Probablity distributions evolving in time
for k in range(0, n):
    plt.plot(x, mod(Psi[k * m])**2, label='$t = {:}$'.format(np.round(t[k*m], 2)))   
    
plt.plot(x, p, 'k--', label="Potential")
plt.legend()
plt.ylim([-0.01, 0.12])
plt.xlabel("X-Position")
plt.title(r"Probablity Distribution Evolution $|\Psi(x)|^2$")
# plt.savefig("fig/FiniteStepProbab.png")
plt.grid()


# Plotting the real and imaginary parts of the wavefunction
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 10))
for k in range(0, 40, 4):
    ax1.plot(x, (Psi[k].real), label='$\Psi(x, t = {:})$'.format(np.round(t[k], 2)))
#     ax1.set_ylim(-0.15, 0.15)
    ax1.set_xlabel('X Position')
    ax1.set_title("Wavefunction Real Part")
    ax2.plot(x, (Psi[k].imag), label='$\Psi(x, t = {:})$'.format(np.round(t[k], 2)))
    ax2.set_xlabel('X Position')
    ax2.set_title("Wavefunction Imaginary Part")
#     ax2.set_ylim(-0.15, 0.15)

ax2.plot(x, p, 'kx--', label="Potential")
ax1.plot(x, p, 'kx--')
ax1.grid()
ax2.grid()
plt.legend()
fig.suptitle("Wavefunction \n Evolution")
# plt.savefig("fig/FiniteStepWaveFun.png")

# Checking normalization of the wavefunction
probab = []
for i in range(len(t)):
    probab.append(np.sum(mod(Psi[i])**2))
    

plt.plot(t[:-1], probab[:-1])  
plt.ylim([0.9, 1.1])
plt.xlabel("Time $t$")
plt.ylabel("Total Probability $\sum_{x=0}^{L} |\Psi(x)|^2$")
plt.title(r"Wavefunction Normalization Check")

plt.grid()

# Running the following cell would render an animation of wavefunction evolution 
# Takes considerable time to run

# animate()