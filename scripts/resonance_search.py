from __future__ import division
from imports import *
import nhqm.bases.mom_space as mom
from nhqm.problems import He5
import nhqm.calculations.QM as calc
from helpers import progress

problem = He5.problem
steps = 50
lowest_energy = sp.empty(steps)
V0s = sp.linspace(-70, 0, steps)
l = 0
j = .5
k_max = 4
order = 20
step_size = k_max / order

plt.title(r"He5, $l = {0}$, $j = {1}$".format(l, j))
plt.xlabel(r"Potential well depth V0 / MeV")
plt.ylabel(r"Ground state energy E / MeV")

for (i, V0) in enumerate(progress(V0s)):
    problem.V0 = V0
    H = calc.hamiltonian(mom.H_element, args=(step_size, problem, l, j), order=order)
    energy, eigvecs = calc.energies(H)
    lowest_energy[i] = energy[0]
    
plt.plot(V0s, lowest_energy)
plt.show()
