from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases import many_body_uncoupled_HO as mb
from nhqm.bases import two_body_interaction as n_n
from nhqm.plot_helpers import *

problem = He5 
order = 10
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 4

n_n.V0 = -1788


Q = osc.QNums(l=1, j=1.5, n=range(order))
H = osc.hamiltonian(order, problem, Q)
eigvals, eigvecs = energies(H)


##################################


Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
             n=range(order),
             m=[-1.5, -0.5, 0.5, 1.5])

mb_H = mb.hamiltonian(Q, eigvals, eigvecs, order, num_particles = 2, verbose=True)
mb_eigvals, mb_eigvecs = energies(mb_H)
print "lowest energy:", mb_eigvals[0]

plt.figure(1)
plt.clf()
plot_poles(mb_eigvals, problem.mass)
#plt.plot(sp.real(mb_eigvals), sp.imag(mb_eigvals), 'ko')

k = sp.zeros((len(points), len(points)),complex)
for i, ks1 in enumerate(points):
    for j, ks2 in enumerate(points):
        k[i,j]=sp.sqrt(ks1 ** 2 + ks2 ** 2)
plt.plot(sp.real(k),sp.imag(k),'og')
        

#plt.show()
