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
contour = calc.gauss_contour(vertices=[0, k_max], order=order)


for (i, V0) in enumerate(progress(V0s)):
    problem.V0 = V0
    args = (problem, l, j)
    H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
    energy, eigvecs = calc.energies(H)
    lowest_energy[i] = energy[0]

plt.title(r"He5, $l = {l}$, $j = {j}$".format(l=l, j=j))
plt.xlabel(r"Potential well depth V0 / MeV")
plt.ylabel(r"Lowest energy E / MeV")  
plt.plot(V0s, lowest_energy)
plt.show()
