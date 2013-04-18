from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour

order = 50
problem = H_atom.problem

Q = osc.QNums(l=0, j=.5, n=range(order))

H = osc.hamiltonian(order, problem, Q)
eigvals, eigvecs = energies(H)
print "Harmonic Oscillator", eigvals[0], problem.units

contour = gauss_contour([0, 5], order)
H = mom.hamiltonian(contour, problem, Q)
eigvals, eigvecs = energies(H)
print "Mom Space", eigvals[0], problem.units