# Computational Physics Tools

This repository is a collection of tools I wrote to learn implementation of some specific computational skills by applying them to some simple (and some not so simple) physical problem. A lot of this comes from refining the scripts I wrote for my 'Computational Physics' and 'Computational Astrophysics and Geophysics' courses. 

Current version has the following: 

1. [Py] Parallelizing Python Code (to solve for Kepler orbits)
2. [Py] Iterative Crank Nicholson Integrator (to solve Schrodinger Equation for arbitrary potentials + generate animations of probability distribution)
3. [Jl] Adaptive Step-Size Runge Kutta (to solve for White-Dwarf Mass Radius Relationship + julia native multithreading)
4. [Py] FFT+Jacobi solver with limited multigrid functionality (to solve for a hydrostatic equilibrium PDE system with small angular momentum)

In the pipeline for near future:

1. [CUDA] GPU Acceleration (to solve N-body problems with C)
2. [Jl] Adaptive Mesh Refinement (to solve, hopefully, some NR system) 
3. [Jl] Metropolis Algorithm (for solving Feynman Path Integral)
4. [Py+Blender] Finite Element Analysis 
