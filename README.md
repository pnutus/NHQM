NHQM – Non-Hermitian Quantum Mechanics
======================================

Library for solving spherically symmetric quantum-mechanical problems, with a focus on using a complex-momentum basis.

Examples
--------
The `scripts` folder contains examples.


Usage
=====

To use the library, choose a problem, a basis to expand the problem in and values for the quantum numbers `l` and `j`. Then you just go:

    problem = He5
    quantum_numbers = QuantumNumbers(l=1, j=1.5)
    basis = MomentumBasis(contour)
    states = solve(problem, quantum_numbers, basis)


Problems
--------

A problem is defined as potential, a characteristic mass and a factor that converts the energy from units of the problem to units of eV.
 
Two problems are included: the Hydrogen atom and the He5 core.
 
Bases
-----

A basis is used to expand in

Two bases expansion methods are included: The isotropic harmonic oscillator and discretized momentum space.
