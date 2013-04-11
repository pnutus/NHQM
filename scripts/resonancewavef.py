from __future__ import division
from imports import *
import nhqm.bases.mom_space as mom
from nhqm.problems import He5
import nhqm.calculations.QM as calc


problem = He5.problem
problem.V0 = -52.3
k_max = 4
order = 20
step_size = k_max / order
l = 1
j = 1.5

H = calc.hamiltonian(mom.H_element, args=(step_size, problem, l, j), order=order)
energy, eigvecs = calc.energies(H)
basis_function = mom.gen_basis_function(step_size, problem, l=l, j=j)
wavefunction = calc.gen_wavefunction(eigvecs[:,0], basis_function)
wavefunction = calc.normalize(wavefunction, 0, 10, weight= lambda r: r**2)

r = sp.linspace(1e-1, 10, 200)

plt.title("Potential, ground energy (and wavefunction, not to scale) \n\
$l = {0}$, $j = {1}$ and $V_0 = {2}$ MeV".format(l, j, problem.V0))
plt.plot(r, l*(l + 1)/ r**2 + problem.potential(r, l, j))
plt.plot(r, .1*r**2 * calc.absq(wavefunction(r)) + energy[0])
plt.plot(r, sp.zeros(len(r)) + energy[0])
plt.axis([0, 10, -0.01, 0.04])
plt.ylabel(r'$r^2|\Psi(r)|^2$')
plt.xlabel(r'$r$ / fm')
plt.show()
