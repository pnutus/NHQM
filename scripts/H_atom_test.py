from __future__ import division
from imports import *
from nhqm.solving import solve
from nhqm.quantum_numbers import QuantumNumbers
from nhqm.bases.momentum import MomentumBasis, gauss_contour
from nhqm.bases.harm_osc import HarmOscBasis
from nhqm.problems import HydrogenAtom
from nhqm.helpers.quantum import absq

# This example calculates the ground state energy of the hydrogen atom
# using by expanding in the harmonic oscillator and momentum bases.
# We also throw in the analytical solutions for good measure.
# The Momentum basis is the worst, as its eigenstates are free particles.

# Setup - try varying these parameters and guess what they do!

problem = HydrogenAtom
basis_state_count = 20
k_max = 7
l = 0
j = 0.5
quantum_numbers = QuantumNumbers(l, j)

# We generate basis objects to use later

contour = gauss_contour([0, k_max], basis_state_count)
mom = MomentumBasis(contour)
osc = HarmOscBasis(basis_state_count, problem.mass, problem.HO_omega)
bases = [osc, mom]

# For plotting

r = sp.linspace(0, 10, 100)

# For each basis, we solve the problem, print the ground state energy
# and plot the probability distribution.

for basis in bases:
    states = solve(problem, quantum_numbers, basis)
    ground_state = states[0]
    
    print basis.name
    print "Ground state energy:", ground_state.energy, problem.energy_units
    
    plt.plot(r, absq(ground_state.wavefunction(r) * r) , label=basis.name )

# Analytical "correct" answer

print "Analytical"
print "Ground state energy:", -0.5, "Hartrees"

def analytical_ground_state_wavefunction(r):
    return 2 * sp.exp(-r)

plt.plot(r, absq(analytical_ground_state_wavefunction(r) * r), 
        label="Analytical Wavefunction")

# A pretty title and legend and we're done!

plt.title("Hydrogen atom ground state radial probability distribution \n"
          "with $l = {1}$, $j = {2}$ "
          "and a ${3} \\times\\,{3}$ matrix"
          .format("", l, j, basis_state_count))
plt.legend()
plt.show()