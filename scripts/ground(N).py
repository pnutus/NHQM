from __future__ import division
from imports import *
import nhqm.bases.mom_space as mom
from nhqm.problems import He5
import nhqm.calculations.QM as calc

problem = He5.problem
problem.V0 = -70
steps = 20
lowest_energy = sp.empty(steps)
orders = sp.arange(steps) + 1
l = 0
j = .5

for (i, order) in enumerate(orders):
    step_size = 2.5 / order
    H = calc.hamiltonian(mom.H_element, \
                args=(problem, l, j, step_size), order=order)
    energy, eigvecs = calc.energies(H)
    lowest_energy[i] = energy[0] # MeV
    print order
    
plt.plot(orders, lowest_energy)
plt.title(r"He5, $l = {0}$, $j = {1}$, $V_0 = {2}$MeV".format(l, j, problem.V0))
plt.xlabel(r"Matrix dimension N ($N \times N$)")
plt.ylabel(r"Ground state energy / MeV")
plt.show()