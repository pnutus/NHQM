NHQM – Non-Hermitian Quantum Mechanics
======================================

Library for calculating energies and wavefunctions for spherically symmetric quantum-mechanical problems. Sometime in the future it might even live up to its name!

Examples
--------
 `example.py` contains examples.


Structure
=========

The library is divided into three parts, problems, methods and hamiltonian calculations. Each separated from the other so that you can choose a problem, which method to use and how to calculate the hamiltonian.

Problems
--------

A problem is defined as potential, a characteristic mass and a factor that converts the energy from units of the problem to units of eV.
 
Two problems are included: the Hydrogen atom and the He5 core.
 
Methods
-------

A method is a way of solving a problem, primarily basis expansion. The basis expansion methods are essentially a function that takes row and column and returns the corresponding matrix element in the hamiltonian.

Two bases expansion methods are included: The isotropic harmonic oscillator and discretized momentum space.

Hamiltonian Calculations
------------------------

A hamiltonian calculation is a way to calculate the hamiltonian. The simplest one is serial -- one at a time.

Will add ways to do this in parallel.