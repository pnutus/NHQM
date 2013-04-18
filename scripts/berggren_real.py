from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour

problem = He5  
order = 30
problem.V0 = -70.
k_max = 4

Q = mom.QNums(l=1, j=1.5, k=range(order))

contour = gauss_contour([0, k_max], order)
H = mom.hamiltonian(contour, problem, Q)
eigvals, eigvecs = energies(H)
print "Berggren on the real line, lowest energy:", eigvals[0]