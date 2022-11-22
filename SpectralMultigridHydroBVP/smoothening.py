import numpy as np
from poissonFFT import *
from myParam import *
import timeit

'''
Smoothening function
Input: Grid Dimension, No of iterations desired,
       optional - Plot Initial State (off by default), Plot final state (off by default), guess value for density (constant by default)
Output: Pressure (complex array), Density (complex array), Potential  (complex array), List of max residue at each iteration
'''
def smooth(N, iterations, initGraph=False,  finGraph=False, guessGiven = 0):
    
    #########################
    ## Initializing Arrays ##
    #########################
    
    x = np.linspace(-(R*1.2), (R*1.2), N)
    y = np.linspace(-(R*1.2), (R*1.2), N)    
    delta = abs(x[-1]-x[0])/N                  # Spatial and Time Step
    rhoGuess = np.zeros((N, N), dtype=complex) # Empty array for storing guess values of density
    
    #########################
    ## Defining Functions ###
    #########################
    
    '''
    Populates the initial guess matrix with density values
    Input: Radial Distance from the center (float)
    Output: Density (Scalar)    
    '''
    def rhoG(r):
        # Defining Boundary
        if abs(r - R) <= delta:
            return rhoB
        elif r < R:
            return rhoB

        # Defining Initial Guess
        else: 
            return 0

    '''
    Computes the residue for one Jacobi Relaxation step
    Input: x-index, y-index, Guess Solution for Pressure (complex matrix of scalars), Guess Solution for potential (complex matrix of scalars)
    Output: Derivative term for Jacobi relaxation (tuple)
    '''
    def residue(i, j, pGuess, phiSol):

        # Coordinate Converstion of term 2
        r = np.sqrt(x[i]**2 + y[j]**2)
        theta = np.arctan2(y[j], x[i])

        # Gradient vector of Pressure
        term_1 = np.array([pGuess[i+1, j] - pGuess[i, j],
                           pGuess[i, j+1] - pGuess[i, j]])/delta

        # Gradient vector of gravitational potential 
        term_2 = np.array([(phiSol[i+1, j] - phiSol[i, j]),
                           (phiSol[i, j+1] - phiSol[i, j])])/delta

        # Gradient vector of centrifugal force    
        term_3 = np.array([2 * r * np.sin(theta)**2,
                           r**2 * np.sin(2*theta)])

        # Component transformation matrix (r, theta) -> (x, y)
        term_3 = pol2cart(term_3)

        return term_1 + (A*term_2 - B*term_3) * rhoGuess[i, j]
    

    '''
    Computes one Relaxation cycle using pressure residue on the grid
    Input: Pressure Guess (complex array), potential guess (complex array)
    Output: Improved Pressure Approximation (complex array), Residue on the grid (complex array) 
    '''
    def JacobiRelaxation(pGuess, phiSol):

        # Initializing Arrays for Pressure
        pStep = np.ones((N, N), dtype=tuple) * pGuess  # Populating first step pressure with guess values
        pSol  = np.zeros((N, N), dtype=tuple)  
        res  = np.zeros((N, N), dtype=tuple)
        
        for i in range(N-1):
            for j in range(N-1):
                r = np.sqrt(x[i]**2 + y[j]**2)
                if abs(r - R) <= delta:
                    pSol[i, j] = pStep[i, j]
                  
                if r < R:
                    res[i,j] = residue(i, j, pStep, phiSol)
                    pSol[i, j] = pStep[i, j] + res[i,j]*delta
                    
        pStep = pSol

        for m in range(N):
            for n in range(N):
                pGuess[m, n] = np.mean(pSol[m,n])
                res[m, n] = np.mean(res[m,n])

        return pGuess, res
    
    '''
    Solve the PDE System iteratively using FFT for Poisson and Jacobi for Pressure
    Input: Relaxation Steps (integer), Guess value for Pressure (complex array), Guess Value for Potential (complex array)
    Output: Relaxed Pressure (complex array), Density (complex array), Potential (complex array), list of maximum residue at each iteration 
    '''
    def varUpdate(steps, pGuess, phiSol):

        maxErList = []
        k = 1 # loop counter
        while k <= steps:

            pGuess, error = JacobiRelaxation(pGuess, phiSol)
            rhoGuess = np.sqrt(pGuess)
            phiSol = PoissonSolveFFT2D(rhoGuess, N, N) 
            
            maxErList.append(np.max(error).real)
            k += 1 # increment

        return pGuess, rhoGuess, phiSol, maxErList
    
    ####################################
    #### End of function definitions ###
    ####################################
    
    # Checking if a density distribution has already been provided
    if np.array_equal(guessGiven, 0):
        # Populating Solution Array with Initial Guess
        for i in range(len(x)):
            for j in range(len(y)):
                rhoGuess[i, j] = rhoG(np.sqrt(x[i]**2 + y[j]**2))
    else:
        rhoGuess = guessGiven
    
    # Converting density to pressure using Equation of state
    pGuess = rhoGuess**2

    # Guess Solution for Phi
    phiSol = PoissonSolveFFT2D(rhoGuess, N, N)
    
    
    # If user commands, the following snippet plots the initial guess values of the variables
    if initGraph == True:
        
        import matplotlib.pyplot as plt
           
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 7))

        ax2.imshow(pGuess.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax2.set_title("Pressure")
        ax1.imshow(rhoGuess.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax1.set_title("Density")
        ax3.imshow(phiSol.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax3.set_title("Potential")
        fig.suptitle('Initial Guess Distribution \n Grid Size = %i' %N, size=20) 
        plt.savefig('fig/preCon.png')
        plt.show()
        
    start_time = timeit.default_timer()  
    # This line computes the solution after relaxation
    p, rho, phi, er = varUpdate(iterations, pGuess, phiSol) 
    elapsed = timeit.default_timer() - start_time  #Computing the runtime of relaxation
        
    # If user commands, the following snippet plots the arrays computed after relaxation
    if initGraph == True:
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 7))

        ax2.imshow(p.real, cmap='inferno', 
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax2.set_title("Pressure")
        ax1.imshow(rho.real, cmap='inferno', 
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax1.set_title("Density")
        ax3.imshow(phi.real, cmap='inferno', 
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax3.set_title("Potential")
        fig.suptitle('Output after %i steps' %iterations, size=20) 
        plt.savefig('fig/postCon.png')
        plt.show()
    
    print(iterations, 'relaxation steps on ', N,' X ', N, ' grid took ', np.round(elapsed, 2), ' sec')
    print('Maximum residue on the grid is ', er[-1])
    
    return p, rho, phi, er

'''
Plots the a heatmap of error on the grid after smoothening has been performed
Input: Error matrix as generated by smoothening function
Output: None - just displays a graph
(This redundant function is defined as a quick fix for some datatype compatibilities)
'''
def plotEr(res):
    from matplotlib.pyplot import colorbar, imshow, title
    N = (len(res))
    resGrid = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            resGrid[i, j] = abs(res[i, j].real)
    a = imshow(resGrid)
    colorbar(a)
    title("Heat Map of Residues on Grid")