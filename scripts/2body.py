from __future__ import division
from imports import *
from nhqm.solution import *
from nhqm.quantum_numbers import QuantumNumbers
from nhqm.bases.mom_space import MomSpaceBasis
from nhqm.bases.harm_osc import HarmOscBasis
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq
from nhqm.bases.contours import gauss_contour

problem = He5
problem.V0 = -70.
basis_state_count = 20
k_max = 7
quantum_numbers = QuantumNumbers(l=0, j=0.5)

contour = gauss_contour([0, k_max], basis_state_count)
mom = MomSpaceBasis(contour)
osc = HarmOscBasis(basis_state_count, problem.mass, problem.HO_omega)
bases = [osc, mom]

for basis in bases:
    print basis.name
    sp_states = solve(problem, quantum_numbers, basis)
    ground_state = sp_states[0]
    print "Ground state energy:", ground_state.energy
    r = sp.linspace(0, 10, 100)
    plt.title(basis.name)
    plt.plot(r, abs(ground_state.wavefunction(r)*r))

plt.show()