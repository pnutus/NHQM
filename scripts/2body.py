from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq, energies
from nhqm.bases.gen_contour import gauss_contour
from collections import namedtuple

problem = He5
problem.V0 = -70.
bases = [osc, mom,]
basis_size = 50
k_max = 7

contour = gauss_contour([0, k_max], basis_size)

QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=0, j=.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))
for basis in bases:
    print basis.name
    eigvals, eigvecs = basis.solution(contour, problem, Q)
    if basis == osc:
        wavef = basis.gen_wavefunction(eigvecs[:,0])
    else:
        wavef = mom.gen_wavefunction(eigvecs[:,0], contour, Q)

    print "Lowest energy:", eigvals[0]
    r = sp.linspace(0, 10, 100)
    plt.title(basis.name)
    plt.plot(r, abs(wavef(r)*r))

plt.show()