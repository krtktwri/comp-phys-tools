from numpy import fft, zeros, cos, sin, pi, array # 
from myParam import * 

'''
Converts polar coordinates to cartesian coordinates (in 2D)
Input: Polar coordinate components (as a numpy array )
Output: Cartesian coordinate components (as a numpy array)
'''
def pol2cart(k):
    x = k[0] * cos(k[1])
    y = k[0] * sin(k[1])

    return array([x, y])

'''
Solves Poisson Equation on a 2 dimensional J x L grid using Fast Fourier Transform 
Input  : Rho (complex array), x-span (int), y-span (int)
Output : Phi array J x L (dtype = complex)
'''

def PoissonSolveFFT2D(rho, J, L):
    
    delta = (2.4 * R)/len(rho)
    rhoHat = fft.fft2(rho)
    phiHat = zeros((J, L), dtype=object)

    for m in range(J):
        for n in range(L):
            denom = 2 * (cos(2 * pi * m/J) + cos(2 * pi * n/L) - 2)
            if denom != 0:
                phiHat[m,n] = (rhoHat[m,n] * delta**2) / denom

    phi = fft.ifft2(phiHat)
    
    return phi

