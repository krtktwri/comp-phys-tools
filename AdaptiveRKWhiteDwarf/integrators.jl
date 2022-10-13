# This code uses CubicSplines.jl (to perform interpolation while inverting the exact equation of state)
# This package do not generally come with Julia installations.

# In Julia REPL, press ] to enter the package mode and type "add CubicSpline" to resolve dependancies
using CubicSplines
include("equations.jl")


# Initializing Arrays
den = range(0.1, 1e12, 2*10^7);     # Array with 10^7 density data points
pre = P_exact.(den);                # Mapping exact EOS to compute pressure for different densities
rho_spline = CubicSpline(pre, den); # Creating a Cubic Spline interpolation using the table populated in the previous line

"
Defining an adaptive step size RK4 integrator to solve 
for mass radius relationship using exact Equation of State's spline interpolation
"
function MassRadius_exact(rho0)
       
    dr = 100.0                       #Spatial step in meters
    e = 1e22

    # Initial Conditions
    r = [1.0];
    rho = [rho0];
    M = [((4/3)π*r[1]^3) * rho0];
    P = [P_exact(rho0)];
    delta = [0.0];
    
    # RK4 Integrator
    
    i = 1
    
    while  P[i] >= 200
               
        rhoLocal = rho_spline(P[i])

        kP1 = dr * dP_dr(r[i], M[i], rhoLocal)
        kM1 = dr * dM_dr(r[i], rhoLocal)
                    
        kP2 = dr * dP_dr(r[i] + dr/2, M[i] + kM1/2, rhoLocal)
        kM2 = dr * dM_dr(r[i] + dr/2, rhoLocal)
    
        kP3 = dr * dP_dr(r[i] + dr/2, M[i] + kM2/2,  rhoLocal)
        kM3 = dr * dM_dr(r[i] + dr/2, rhoLocal)

        kP4 = dr * dP_dr(r[i] + dr, M[i] + kM3, rhoLocal)
        kM4 = dr * dM_dr(r[i] + dr, rhoLocal)
        
        append!(delta, abs(kM4- kM1))

        if delta[i] <= e
            append!(r, r[i] + dr)
            append!(P, P[i] + (kP1 + 2kP2 + 2kP3 + kP4)/6)
            append!(M, M[i] + (kM1 + 2kM2 + 2kM3 + kM4)/6)
        
        else

            rhoLocal = rho_spline(P[i])            
            dr = dr * (e/delta[i])^(1/5)
            
            kP1 = dr * dP_dr(r[i], M[i], rhoLocal)
            kM1 = dr * dM_dr(r[i], rhoLocal)
                        
            kP2 = dr * dP_dr(r[i] + dr/2, M[i] + kM1/2, rhoLocal)
            kM2 = dr * dM_dr(r[i] + dr/2, rhoLocal)
        
            kP3 = dr * dP_dr(r[i] + dr/2, M[i] + kM2/2,  rhoLocal)
            kM3 = dr * dM_dr(r[i] + dr/2, rhoLocal)
    
            kP4 = dr * dP_dr(r[i] + dr, M[i] + kM3, rhoLocal)
            kM4 = dr * dM_dr(r[i] + dr, rhoLocal)
            
            append!(r, r[i] + dr)
            append!(P, P[i] + (kP1 + 2kP2 + 2kP3 + kP4)/6)
            append!(M, M[i] + (kM1 + 2kM2 + 2kM3 + kM4)/6)
        
        end
        
        i += 1
    
    end
    
    return M[end-1], r[end-1]

end

"
Defining an adaptive step size RK4 integrator to solve 
for mass radius relationship using approximate Equation of State 
"

function MassRadius_approx(rho0)
       
    dr = 100.0                       #Spatial step in meters
    e = 1e22

    # Initial Conditions
    r = [1.0];
    rho = [rho0];
    M = [((4/3)π*r[1]^3) * rho0];
    P = [P_exact(rho0)];
    delta = [0.0];
    
    # RK4 Integrator
    
    i = 1
    
    while  P[i] >= 200
               
        rhoLocal = rho_approx(P[i])

        kP1 = dr * dP_dr(r[i], M[i], rhoLocal)
        kM1 = dr * dM_dr(r[i], rhoLocal)
                    
        kP2 = dr * dP_dr(r[i] + dr/2, M[i] + kM1/2, rhoLocal)
        kM2 = dr * dM_dr(r[i] + dr/2, rhoLocal)
    
        kP3 = dr * dP_dr(r[i] + dr/2, M[i] + kM2/2,  rhoLocal)
        kM3 = dr * dM_dr(r[i] + dr/2, rhoLocal)

        kP4 = dr * dP_dr(r[i] + dr, M[i] + kM3, rhoLocal)
        kM4 = dr * dM_dr(r[i] + dr, rhoLocal)
        
        append!(delta, abs(kM4- kM1))

        if delta[i] <= e
            append!(r, r[i] + dr)
            append!(P, P[i] + (kP1 + 2kP2 + 2kP3 + kP4)/6)
            append!(M, M[i] + (kM1 + 2kM2 + 2kM3 + kM4)/6)
        
        else
            
            dr = dr * (e/delta[i])^(1/5)
            
            kP1 = dr * dP_dr(r[i], M[i], rhoLocal)
            kM1 = dr * dM_dr(r[i], rhoLocal)
                        
            kP2 = dr * dP_dr(r[i] + dr/2, M[i] + kM1/2, rhoLocal)
            kM2 = dr * dM_dr(r[i] + dr/2, rhoLocal)
        
            kP3 = dr * dP_dr(r[i] + dr/2, M[i] + kM2/2,  rhoLocal)
            kM3 = dr * dM_dr(r[i] + dr/2, rhoLocal)
    
            kP4 = dr * dP_dr(r[i] + dr, M[i] + kM3, rhoLocal)
            kM4 = dr * dM_dr(r[i] + dr, rhoLocal)
            
            append!(r, r[i] + dr)
            append!(P, P[i] + (kP1 + 2kP2 + 2kP3 + kP4)/6)
            append!(M, M[i] + (kM1 + 2kM2 + 2kM3 + kM4)/6)
        
        end
        
        i += 1
    
    end
    
    return M[end-1], r[end-1]

end
