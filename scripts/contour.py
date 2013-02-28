from __future__ import division
from imports import *
from plotting import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc

problem = He5.problem   
order = 100
l = 1
j = 1.5
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.05
k_max = 2.5

contour = mom.gen_simple_contour(peak_x, peak_y, k_max, order)

args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)
print "Berggren lowest energy:", eigvals[0]
basis_function = mom.gen_basis_function(args)
wavefunction = calc.gen_wavefunction(eigvecs[:,0], basis_function, contour)
r = sp.linspace(0, 10)
plt.figure(0)
plt.plot(r, r**2 * calc.absq(wavefunction(r)))
plt.figure(1)
plot_poles(eigvals, problem.mass)
plt.figure(2)
plot_wavefunctions(contour, eigvecs)
plt.show()