from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.bases import berggren as berg
from nhqm.problems import He5
from nhqm.calculations import QM as calc

problem = He5.problem   
order = 30
l = 1
j = 1.5
problem.V0 = -70.

contour = sp.linspace(0, 2.5, order + 1)
args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)
print "MomSpace:", eigvals[0]
H = calc.contour_hamiltonian(berg.H_element, contour, args)
eigvals, eigvecs = calc.energies(H)
print "Berggren:", eigvals[0]