from __future__ import division
from imports import *
from nhqm.helpers.plotting import *
from nhqm.bases.momentum import MomentumBasis, triangle_contour
from nhqm.problems import Helium5
from nhqm.solve import solve
from nhqm.quantum_numbers import QuantumNumbers
from nhqm.helpers.timing import progress

# This example plots the resonance pole as the potential well
# gets more and more shallow. The pole is first a bound state on
# the imaginary axis, then gets less bound and finally moves into the 
# fourth quadrant and becomes a resonance.

# Setup - try experimenting with the numbers!

quantum_numbers = QuantumNumbers(l=1, j=1.5)
basis_state_count = 20
triangle_state_count = 16
peak_x = 0.17
peak_y = -0.07
k_max = 2.5
contour = triangle_contour(peak_x, peak_y, k_max, triangle_state_count, 
                           basis_state_count - triangle_state_count)
momentum_basis = MomentumBasis(contour)

# Parameters for the potential strength variation

start_V0 = -55 #MeV
end_V0   = -47 #Mev
steps = 20
V0s = sp.linspace(start_V0, end_V0, steps)
ks = sp.empty(steps, complex)


def positive_if_imaginary(k):
    """
    Hack to fix rounding errors that make bound states appear
    on the negative imaginary axis.
    """
    if sp.real(k) < 1e-4 and abs(sp.imag(k)) > 0.02:
        return 1j*abs(k)
    else:
        return k

# Solve the problem and save the momentum (k) values for 
# every resonance state

for (i, V0) in enumerate(progress(V0s)):
    problem = Helium5(V0 = V0) 
    states = solve(problem, quantum_numbers, momentum_basis)
    resonance = find_resonance_state(states)
    k = sp.sqrt(2 * problem.mass * resonance.energy)
    ks[i] = positive_if_imaginary(k)

# Plot the resonance ks together with the contour

plot_contour(contour)
plt.plot(sp.real(ks), sp.imag(ks), 'ro')
plt.title("Movement of the resonance pole as the potential s")
plt.show()
