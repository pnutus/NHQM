NHQM – Non-Hermitian Quantum Mechanics
======================================

Library for calculating energies and wavefunctions for spherically symmetric quantum-mechanical problems. Sometime in the future it might even live up to its name!

Examples
--------
The `scripts` folder contains examples.


Structure
=========

The library is divided into three parts, problems, bases and calculations. Each separated from the other so that you can choose a problem, which method to use and how to calculate the hamiltonian.

Problems
--------

A problem is defined as potential, a characteristic mass and a factor that converts the energy from units of the problem to units of eV.
 
Two problems are included: the Hydrogen atom and the He5 core.
 
Bases
-----

A basis module represents a basis to expand in. The basis expansion methods are essentially a function that takes row and column and returns the corresponding matrix element in the hamiltonian.

Two bases expansion methods are included: The isotropic harmonic oscillator and discretized momentum space.

Calculations
------------

Calculations are ways to calculate the hamiltonian, wavefunctions, etc. The simplest one is serial -- one at a time.

Will add ways to do this in parallel.