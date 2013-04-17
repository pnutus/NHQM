from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import triangle_contour
from nhqm.bases import many_body as mb
from nhqm.bases import two_body_interaction as n_n
from nhqm.plot_helpers import *

problem = He5.problem   
order = 3*3
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 4

n_n.V0 = -100

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

contour = triangle_contour(peak_x, peak_y, k_max, order/3)

Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
             m=[-1.5, -0.5, 0.5, 1.5], 
             k=range(len(contour[0])))
             
H = mom.hamiltonian(contour, problem, Q)
eigvals, eigvecs = energies(H)

mb_H = mb.hamiltonian(Q, eigvals, eigvecs, contour, num_particles = 2, verbose=True)
mb_eigvals, mb_eigvecs = energies(mb_H)
print "lowest energy:", eigvals[0]
plt.plot(sp.real(eigvals), sp.imag(eigvals), 'k')
plt.show()