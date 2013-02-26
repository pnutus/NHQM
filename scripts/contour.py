from __future__ import division
from imports import *
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

def plot_poles(energies, mass):
    ks = sp.sqrt(2*mass*energies)
    plt.plot(sp.real(ks), sp.imag(ks), 'o')
    plt.show()



contour = mom.gen_simple_contour(peak_x, peak_y, k_max, order)

args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)
print "Berggren lowest energy:", eigvals[0]
plot_poles(eigvals, problem.mass)