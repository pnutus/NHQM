NHQM – Non-Hermitian Quantum Mechanics
======================================

Library for solving spherically symmetric quantum-mechanical problems through
basis expansion, with a focus on a complex-momentum basis.

Examples
--------
The `scripts` folder contains examples. Be sure to look at those. 


Usage
=====

To use the library, choose a problem, a basis to expand the problem in and values for the quantum numbers `l` and `j`. Import what you need:

    from nhqm.problems import Helium5
    from nhqm.quantum_numbers import QuantumNumbers
    from nhqm.bases.momentum import MomentumBasis, gauss_contour
    from nhqm.solve import solve

Then you just:

    problem = Helium5()
    quantum_numbers = QuantumNumbers(l=1, j=1.5)
    basis = MomentumBasis(contour)
    states = solve(problem, quantum_numbers, basis)

The states obtained are solutions to the time-independent Schrödinger equation,
i.e. stationary states. These can be used as a basis for a many body problem.

Problems
--------

A problem is an object with a potential and a characteristic mass. It has some other parameters as well for convenience.
 
Two problems are included: the Hydrogen atom and the He5 core.
 
Bases
-----

To solve a problem, we expand it in a basis. The basis is an object with methods that can create a Hamiltonian matrix and generate wavefunctions from the eigenvectors of the Hamiltonian.

Two bases are included: The isotropic harmonic oscillator and the momentum basis (free particle basis). The momentum basis can use complex momenta by integrating along a complex contour, yielding a _non-hermitian_ (but symmetric) hamiltonian with possibly complex eigenvalues. The eigensolutions to such a Hamiltonian form the so-called _Berggren basis_ which can be used as a single-particle basis when working with many-particle systems.
