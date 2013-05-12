from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases import many_body_uncoupled_HO as mb
from nhqm.bases import two_body_interaction_HO as n_n
from nhqm.bases.gen_contour import triangle_contour
from nhqm.plot_helpers import *

from nhqm.bases import mom_space as mom

problem = He5 
order = 50
problem.V0 = -47
n_n.V0 = -180

Q = osc.QNums(l=1, j=1.5, n=range(order))
H = osc.hamiltonian(order, problem, Q)
eigvals, eigvecs = energies(H)

#print eigvals

##################################


Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
             E=range(order),
             m=[-1.5, -0.5, 0.5, 1.5])

mb_H = mb.hamiltonian(Q, eigvals, eigvecs, num_particles = 2, verbose=True)
mb_eigvals, mb_eigvecs = energies(mb_H)
print "lowest energy:", mb_eigvals[0]

#print mb_H