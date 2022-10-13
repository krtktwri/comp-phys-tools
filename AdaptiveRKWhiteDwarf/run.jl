# Calling the integrator scripts. Some overhead due to spline interpolation calculation
include("integrators.jl")

rhoRange = [10^i for i in 7:12];           # Defining central densities in Kg/m^3
const final_index = length(rhoRange)

# Initializing empty arryays
sol_exact = Array{Tuple}(undef,  final_index) 
sol_approx = Array{Tuple}(undef, final_index) 

# Populating with solutions
Threads.@threads for n in 1:final_index
    sol_exact[n] = MassRadius_exact(rhoRange[n])            # Computing solutions with exact Equation of State's spline interpolation
    sol_approx[n] = MassRadius_approx(rhoRange[n])          # Computing solutions with approximate non-relativistic Equation of State
end

# Plotting using PyPlot
# This code uses PyPlot.jl (to perform interpolation while inverting the exact equation of state)
# This package do not generally come with Julia installations.

# In Julia REPL, press ] to enter the package mode and type "add PyPlot" to resolve dependancies
using PyPlot

# Julia's native Plots.jl has a massive compute overhead (for optimization purposes)
# But it reduces the time-to-first-plot, which reduces performance in our small plotting code

# Reshaping the solutions arrays appropriately
M_e, R_e = [x[1] for x in sol_exact], [x[2] for x in sol_exact]
M_a, R_a = [x[1] for x in sol_approx], [x[2] for x in sol_approx]


plot(M_e/1.988e30, R_e/1000, "o--", label="Exact EOS") 
plot(M_a/1.988e30, R_a/1000, "x--", label="Approximate EOS")
xlabel("Mass (in Solar Mass)")
ylabel("Radius (in Kms)")
title("White Dwarf Mass Radius Relationship")
vlines(1.42, 0, 3e4, linestyles="dotted", label="Chandrashekhar Limit")
xlim([0, 2])
legend() 
grid(alpha=0.1)
savefig("MassRadiusRelations.png")