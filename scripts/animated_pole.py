from __future__ import division
from imports import *
from nhqm.plot_helpers import *
from nhqm.bases.momentum import MomentumBasis, triangle_contour
from nhqm.problems import He5
from nhqm.solution import solve
from nhqm.quantum_numbers import QuantumNumbers

# This example plots the resonance pole as the potential well
# gets more and more shallow. The pole is first a bound state on
# the imaginary axis, then gets less bound and finally moves into the 
# fourth quadrant and becomes a resonance.

# Setup - try experimenting with the numbers!

problem = He5  
quantum_numbers = QuantumNumbers(l=1, j=1.5)
basis_state_count = 30
peak_x = 0.17
peak_y = 0.07
k_max = 2.5
contour = triangle_contour(peak_x - 1j*peak_y, k_max, basis_state_count/3)
momentum_basis = MomentumBasis(contour)

# Parameters for the potential strength variation

startV0 = -55
endV0 = -45
steps = 20
V0s = sp.linspace(startV0, endV0, steps)
ks = sp.empty(steps, complex)

# Solve the problem and save the momentum (k) values for 
# every resonance state

def force_positive_imaginary_part(k):
    """
    Hack to fix rounding errors.
    """
    if sp.real(k) < 1e-4 and abs(sp.imag(k)) > 0.02:
        return 1j*abs(k)
    else:
        return k

for (i, V0) in enumerate(V0s):
    problem.V0 = V0
    states = solve(problem, quantum_numbers, momentum_basis)
    resonance = find_resonance_state(states)
    k = sp.sqrt(2*problem.mass*resonance.energy)
    ks[i] = force_positive_imaginary_part(k)

# Plot the momentum poles together with the contour

plot_contour(contour)
plt.plot(sp.real(ks), sp.imag(ks), 'ro')
plt.show()
