import numpy as np                  # For dealing with arrays effectively
from multiprocessing import Pool    # To parallelize the RK4 and Euler solvers 
from functools import partial       # To pool functions with multiple arguments

'''
Computes acceleration for a mass in one-body kepler problem
------
Parameters:
Position (3D Vector)

Returns: 
Acceleration (3D Vector)
'''

def centralForce(r):
    x, y, z = r
    
    r_cube = (x**2 + y**2 + z**2)**1.5 # Norm cubed
    
    ax = -GM*x/r_cube
    ay = -GM*y/r_cube
    az = -GM*z/r_cube
    
    return np.array([ax, ay, az])

'''
Calculates Energy per unit mass
------
Parameters:
Position (3D Vector), Velocity (3D Vector)

Returns: 
Total Energy (per Unit Mass), Kinetic Energy (per Unit Mass), Potential Energy (per Unit Mass)
'''

def scaledE(r, v):
    
    T = 0.5*np.linalg.norm(v)**2   # Kinetic Energy
    V = -GM/np.linalg.norm(r)      # Potential Energy
    
    return T+V, T, V

'''
Calculates Angular Momentum per unit mass
------
Parameters:
Position (3D Vector), Velocity (3D Vector)

Returns: 
Angular Momentum (3D Vector)
'''

def scaledL(r, v):
    
    x, y, z = r      # Unpacking Position
    vx, vy, vz = v   # Unpacking Velocity
    
    # Cross product computation for Angular momentum
    lx, ly, lz = (y*vz - z*vy), (z*vx - x*vz), (x*vy - y*vx)
    
    return np.array([lx, ly, lz])


########################
#### Integrators #######
########################

'''
Runge-Kutta (Fourth Order) Integrator for Second Order ODE
------
Parameters:
Initial Position (Vector), Initial Velocity (Vector), Number of orbits, time-step = 1e5 seconds

Returns: 
Position Vector, Velocity Vector, Time Array
'''

def trajectoryRK4(initPos, initVel, orbits, dt = 1e5):
    
    tf = orbits * period(a, M)   # Final time value
    t = np.arange(0, tf, dt)     # Time array
    N = len(t)                   # Number of integration steps
    
    r = np.zeros((N, 3), float)  # Empty array for position vectors
    v = np.zeros((N, 3), float)  # Empty array for velocity vectors

    r[0] = np.array([initPos])   # Setting initial position from input parameter
    v[0] = np.array([initVel])   # Setting initial velocity from input parameter 
    
    # Runge Kutta Fourth Order Integrator
    for i in range(1, N):

        kv1 = dt * centralForce(r[i-1])
        kr1 = dt * v[i-1] 
        
        kv2 = dt * centralForce(r[i-1] + (kr1/2))
        kr2 = dt * (v[i-1] + (kv1/2))

        kv3 = dt * centralForce(r[i-1] + (kr2/2))
        kr3 = dt * (v[i-1] + (kv2/2))

        kv4 = dt * centralForce(r[i-1] + (kr3)) 
        kr4 = dt * (v[i-1] + kv3)

        #Advancing using the weighted average of four intermediate steps 
        v[i] = v[i-1] + (kv1 + 2*kv2 + 2*kv3 + kv4)/6
        r[i] = r[i-1] + (kr1 + 2*kr2 + 2*kr3 + kr4)/6 
    
    print("Runge-Kutta (4th Order): Finished integrating the trajectory for ", orbits, " orbits")
    return r, v, t

'''
Euler-Cromer Integrator for Second Order ODE
------
Parameters:
Initial Position (Vector), Initial Velocity (Vector), Number of orbits, time-step = 1e5 seconds

Returns:
Position Vector, Velocity Vector, Time Array
'''

def trajectoryEuler(initPos, initVel, orbits, dt = 1e5):
    
    tf = orbits * period(a, M)   # Final time value
    t = np.arange(0, tf, dt)     # Time array
    N = len(t)                   # Number of integration steps
    
    r = np.zeros((N, 3), float)  # Empty array for position vectors
    v = np.zeros((N, 3), float)  # Empty array for velocity vectors

    r[0] = np.array([initPos])   # Setting initial position from input parameter
    v[0] = np.array([initVel])   # Setting initial velocity from input parameter 
    
    # Euler Cromer Integrator
    for i in range(1, N):

        r[i] = r[i-1] + v[i-1]*dt 
        v[i] = v[i-1] + centralForce(r[i])*dt

    
    print("Euler-Cromer: Finished integrating the trajectory for ", orbits, " orbits")   
        
    return r, v, t

# Computes period (in seconds) for Kepler Orbits
def period(r, M):            
    return np.sqrt(4*(np.pi**2) * ((r/a)**3)/GM) 

########################
##### Constants ########
########################
a = 1.521e11
M = 1.989e30                 # in KGs
G = 6.674e-11                # in m^3/kg.s^2
GM = G*M/(a**3)              # Rescaling GM in AUs
yr = period(a, M)            # 1 year in seconds 

########################
## Initial conditions ##
########################

initPos = 1.521e11/a, 0, 0   # in AU
initVel = 0, 29290/a, 0      # in AU/s

if __name__ == '__main__':
   
    # Pooling RK4 and Euler solvers to run both parallely 
    with Pool(3) as p:
    
    #######################
    ###### Solver #########
    #######################
  
      compute1 = p.map(partial(trajectoryRK4, initPos, initVel), [100])
      compute2 = p.map(partial(trajectoryEuler, initPos, initVel), [100])

    #######################
    ##### Unpacking Data ##
    #######################
    
    rRK, vRK, t = compute1[0][0], compute1[0][1], compute1[0][2] 
    rEC, vEC, t = compute2[0][0], compute2[0][1], compute2[0][2]

    '''
    The following commands dump data arrays in respective folders. 
    This is of use if being run on a cluster. 
    One can generate output data, download it using SFTP and proceed with plotting on the local device. 
    '''
    import os
    import shutil

    # os library has issues overwriting existing directories. 
    # This will check if the directory already exists and delte it.
    dir="RKData"
    if os.path.exists(dir):
        shutil.rmtree(dir)

    os.mkdir('RKData')                                # Making a folder for storing RK4 data
    np.savetxt("RKData/rRK.csv", rRK, delimiter=",")
    np.savetxt("RKData/vRK.csv", vRK, delimiter=",")
    np.savetxt("RKData/t.csv", t, delimiter=",")

    
    # os library has issues overwriting existing directories. 
    # This will check if the directory already exists and delte it.
    dir="EulerData"
    if os.path.exists(dir):
        shutil.rmtree(dir)

    os.mkdir('EulerData')                            # Making a folder for storing Euler-Crommer data
    np.savetxt("EulerData/rEC.csv", rRK, delimiter=",")
    np.savetxt("EulerData/vEC.csv", vRK, delimiter=",")
    np.savetxt("EulerData/t.csv", t, delimiter=",")


    #########################
    ### Unpacking Vectors ###
    #########################
    xRK, yRK, zRK = rRK[:, 0], rRK[:, 1], rRK[:, 2]
    xEC, yEC, zEC = rEC[:, 0], rEC[:, 1], rEC[:, 2]

    '''
    Performing Energy and Angular Momentum Checks
    '''
    ERK0 = scaledE(rRK[0], vRK[0])[0]
    ERKf = scaledE(rRK[-1], vRK[-1])[0]

    LRK0 = np.linalg.norm(scaledL(rRK[0], vRK[0]))
    LRKf = np.linalg.norm(scaledL(rRK[-1], vRK[-1]))

    EEC0 = scaledE(rEC[0], vEC[0])[0]
    EECf = scaledE(rEC[-1], vEC[-1])[0]

    LEC0 = np.linalg.norm(scaledL(rEC[0], vEC[0]))
    LECf = np.linalg.norm(scaledL(rEC[-1], vEC[-1]))


    print("RK4 (100 orbits): \n Energy Variation = ", (np.abs((ERK0 - ERKf)/ERK0) * 100), "% \n Ratio of Initial and Final Angular Momentum = ", LRK0/LRKf)
    print("Euler-Cromer (100 orbits): \n Energy Variation = ", (np.abs((EEC0 - EECf)/EEC0) * 100),'% \n Ratio of Initial and Final Angular Momentum = ', LEC0/LECf)

    # A strange observation, when I ran the code, 
    # was that even though the energy is better conserved in RK4 integration, 
    # for some reason Angular momentum is better conserved in Euler Crommer Integration 


    '''
    The following commands are for plotting purposes. 
    They have been commented out in order to make sure the code runs without error in all environments 
    Some people might not have matplotlib or the DLLs might not be correct. Or, some people might choose to run the script on a cluster.
    '''

#     import matplotlib.pyplot as plt
#     import warnings
#     warnings.filterwarnings( "ignore", module = "matplotlib\..*" )

    
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 10))
#     plt.suptitle("Simulation Duration = 100 Orbital Periods")
#     ax1.plot(xRK, yRK, label="RK Trajectory")
#     ax1.set_xlabel('X Position [AU]')
#     ax1.set_ylabel('Y Position [AU]')
#     ax1.set(aspect=1)
#     ax1.set_title("Trajectory Integrated by RK4")


#     ax2.plot(xEC, yEC, label="Euler Trajectory")
#     ax2.set_xlabel('X Position [AU]')
#     ax2.set_ylabel('Y Position [AU]')
#     ax2.set(aspect=1)
#     ax2.set_title("Trajectory Integrated by Euler-Cromer")

