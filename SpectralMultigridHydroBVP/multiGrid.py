from numpy import ones, sqrt, linspace
from myParam import R, rhoB
from smoothening import smooth

'''
Restriction Operator (uses 9 target points)
Input: Numpy Array of size 2N x 2N
Output: Numpy Array scaled down to N x N, new dimension of the matrix
'''
def restriction(test_c):
    
    N = len(test_c)
    newN = int(N/2)
    out = ones((newN, newN), dtype=complex)
    
    for i in range(1, newN-1):
        for j in range(1, newN-1):
            out[i , j] = (1/16) * (4 * test_c[2*i, 2*j] + 
                                      2 * (test_c[2*i, 2*j-1] + test_c[2*i, 2*j+1] + test_c[2*i+1, 2*j] + test_c[2*i-1, 2*j]) 
                                      + (test_c[2*i-1, 2*j-1] + test_c[2*i+1, 2*j-1] + test_c[2*i-1, 2*j+1] + test_c[2*i+1, 2*j+1]))
    
    return out*10, newN

'''
Prolongation Operator (uses 9 target points)
Input: Numpy Array of size N x N
Output: Numpy Array scaled down to 2N x 2N, new dimension of the matrix
'''
def prolongation(test_f):
        
    N = len(test_f)
    newN = int(2*N)
    out = ones((newN, newN), dtype=complex)
    
    for i in range(1, N-1):
        for j in range(1, N-1):
            out[2*i, 2*j]  = test_f[i, j]
            out[2*i+1, 2*j] = 0.5*(test_f[i+1, j] + test_f[i, j])
            out[2*i, 2*j+1] = 0.5*(test_f[i, j+1] + test_f[i, j])
            out[2*i+1, 2*j+1] = 0.25*(test_f[i, j] + test_f[i+1, j] + test_f[i, j+1] + test_f[i+1, j+1])
                        
    return out/10, newN

'''
Prolongation Operator (uses 9 target points) but artificially imposes boundary conditions again after rescaling
Input: Numpy Array of size N x N
Output: Numpy Array scaled down to 2N x 2N, new dimension of the matrix
'''
def prolongImposeBC(coarse):
    fine, newN = prolongation(coarse)
    x = linspace(-(R*1.2), (R*1.2), newN)
    y = linspace(-(R*1.2), (R*1.2), newN)
    delta = abs(x[-1]-x[0])/newN
    
    for i in range(newN):
        for j in range(newN):
            r = sqrt((x[i]**2 + y[j]**2))
            if abs(r - R) <= delta:
                fine[i, j] = rhoB
            elif r > R:
                fine[i, j] = 0
    
    return fine, newN

'''
Restriction Operator (uses 9 target points) but artificially imposes boundary conditions again after rescaling
Input: Numpy Array of size 2N x 2N
Output: Numpy Array scaled down to N x N, new dimension of the matrix
'''
def restrictImposeBC(fine):
    coarse, newN = restriction(fine)
    x = linspace(-(R*1.2), (R*1.2), newN)
    y = linspace(-(R*1.2), (R*1.2), newN)
    delta = abs(x[-1]-x[0])/newN
    
    for i in range(newN):
        for j in range(newN):
            r = sqrt((x[i]**2 + y[j]**2))
            if abs(r - R) <= delta:
                coarse[i, j] = rhoB
            elif r > R:
                coarse[i, j] = 0
    
    return fine, newN

'''
Solves using FFT+Jacobi but with an initial guess arrived at by solving on a coarse grid and prolongating the solution to a finer grid
Input: Desired grid size N, number of relaxation steps, 
       optional - Boolean, Save Graphs (False by default)
'''
def multiGridRelax(N, steps, graphSave = False):
    
    N_guess = int(N/2)
    print("First solving on grid of size N = ", N_guess)
    pSolG, rhoSolG, phiSolG, er = smooth(N_guess, int(0.5*steps), False, False)
        
    print("Now prolongating to desired resolution")
    fine, fineN = prolongImposeBC(rhoSolG)
    pSol, rhoSol, phiSol, er = smooth(fineN, steps, False, False, fine)
    
    if graphSave == True:
    
        import matplotlib.pyplot as plt
           
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

        ax1.imshow(rhoSolG.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax1.set_title("Solution on Coarse Grid N = %i" %N_guess)
        
        ax2.imshow(fine.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax2.set_title("Grid Prolongation to N = %i" %fineN)

        fig1.suptitle('Utilizing Multigrid Machinery', size=20) 
        plt.savefig('fig/preCon_multiGrid.png')
        plt.show()
        
        fig2, (ax4, ax5, ax6) = plt.subplots(1, 3, figsize=(20, 7))
        
        ax4.imshow(rhoSol.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax4.set_title("Density")
        ax5.imshow(pSol.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax5.set_title("Pressure")
        ax6.imshow(phiSol.real, cmap='gnuplot',
                   extent =[-1.2*R, 1.2*R,-1.2*R, 1.2*R])
        ax6.set_title("Potential")
        fig2.suptitle('Final solution on finer grid \n Grid Size = %i' %fineN, size=20) 
        plt.savefig('fig/postCon_multiGrid.png')
        plt.show()
        
    return pSol, rhoSol, phiSol, er   
    