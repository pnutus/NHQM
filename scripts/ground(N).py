from __future__ import division
from imports import *
import basis_mom_space as mom
import problem_He5 as He5
import calculate_serial as calc
from wavefunction import gen_wavefunction
from normalize import normalize, norm

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