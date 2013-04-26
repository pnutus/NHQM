from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour

order = 50
problems = H_atom, He5
mom.integration_order =30


mom.integration_order = 20
mom.integration_range = 10
osc.integration_order = 60

Q = osc.QNums(l=0, j=.5, n=range(order))

for problem in problems:
    print problem.name + ":"
    H = osc.hamiltonian(order, problem, Q)
    eigvals, eigvecs = energies(H)
    print osc.name, eigvals[0], problem.units

    contour = gauss_contour([0, 5], order)
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    print mom.name, eigvals[0], problem.units