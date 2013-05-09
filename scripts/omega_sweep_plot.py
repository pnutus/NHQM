from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.problems import He5, H_atom
from nhqm.calculations import QM as calc

problem = He5
bases = osc
order = 30
k_max = 5
step_size = k_max / order
l = 0
j = 0.5
r = sp.linspace(0, 10, 200)

Q = osc.QNums(l=0, j=.5, n=range(order))

def absq(x):
    return x*sp.conjugate(x)
    
plt.suptitle("Radial probability distribution "
          "with $l = {1}$, $j = {2}$ "
          "and a ${3} \\times\\,{3}$ matrix"
          .format("", l, j, order))


plt.title(problem.name)
plt.ylabel(r'$r^2|\Psi|^2$')
plt.xlabel(r'$r$')
plt.show(block=False)

omega = 1
H = osc.hamiltonian(order, problem, Q)

energy, eigvecs = calc.energies(H)
lowest_energy = energy[0] * problem.eV_factor

plt.legend()
plt.draw()

plt.show() # pause when done