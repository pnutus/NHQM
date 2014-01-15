# NHQM
## Non-Hermitian Quantum Mechanics
Library for solving spherically symmetric quantum-mechanical problems through
basis expansion, with a focus on a complex-momentum basis.

This code was part of a Bachelor thesis project. For more information, see the [thesis](http://publications.lib.chalmers.se/records/fulltext/179709/179709.pdf). Some functionality has been removed during cleanup since then, so a [legacy branch](https://github.com/pnutus/NHQM/tree/legacy) was created for historical purposes.

## Usage

To use the library, choose a problem, a basis to expand the problem in and values for the quantum numbers `l` and `j`. Import what you need:

    from nhqm.solve import solve
    from nhqm.quantum_numbers import QuantumNumbers
    from nhqm.bases.momentum import MomentumBasis, gauss_contour
    from nhqm.problems import Helium5

Then do something like this:

    problem = Helium5()
    quantum_numbers = QuantumNumbers(l=1, j=1.5)
    contour = gauss_contour([0, 10], points_per_segment=10)
    basis = MomentumBasis(contour)
    states = solve(problem, quantum_numbers, basis)

The states obtained are solutions to the time-independent Schrödinger equation, i.e. stationary states. These can be used as a basis for a many body problem. The states contain information you can use:

    ground_state = states[0]
    print ground_state.energy
    plot(ground_state.wavefunction) #pseudo-code

### Examples
The `scripts` folder contains examples. Be sure to run and look at those. 

##Problems

A problem is an object with a potential and a characteristic mass. It has some other parameters as well for convenience.
 
Two problems are included: the Hydrogen atom and the He5 core.
 
##Bases

To solve a problem, we expand it in a basis. The basis is an object with methods that can create a Hamiltonian matrix and generate wavefunctions from the eigenvectors of the Hamiltonian.

Two bases are included: The isotropic harmonic oscillator and the momentum basis (free particle basis).

### The Complex-Momentum Basis

The momentum basis can use complex momenta by integrating along a complex contour, yielding a _non-hermitian_ (but symmetric) hamiltonian with possibly complex eigenvalues. The eigensolutions to such a Hamiltonian form the so-called _Berggren basis_ which can be used as a single-particle basis when working with many-particle systems.

## States

When solving a single-particle problem, the result is a list of `SingleParticleStationaryState` objects. These objects each represent a solution to the time-independent Schrödinger equation and contain information about problem, basis, quantum numbers, energy and wavefunction. This object representation can be useful when using the states as a single-particle basis for a many-body system.

## Quantum Numbers

The Quantum Number object is very simple, only containing `l` and `j`.

## Helpers

There are a number of helper functions included, ranging from very important to unimportant. Most important are the `matrix` and `quantum` modules.
