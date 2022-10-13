# Constants
const solarMass = 1.988e30         #in Kilogram
const G = 6.674e-11 

# Differential Equations
dP_dr(r, M, rho) = -(G*M/r^2)*rho;
dM_dr(r, rho) = 4π * r^2 * rho;

# Equations of State
rho_approx(P) = ((P/(3.16e21))^(3/5))*10.0^9;       #Approximate EOS (ignoring special relativity)

function P_exact(rho)          
    x = 0.80 * (rho*1e-9)^(1/3)
    return 1.8e22 * (x * (1+x^2)^(1/2) * (((2*(x^2))/3) - 1) + log(x + (1 + x^2)^(1/2)))
end;
