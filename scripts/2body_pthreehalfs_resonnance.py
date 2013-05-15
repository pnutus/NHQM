from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq, energies
from nhqm.bases.gen_contour import gauss_contour, triangle_contour
from collections import namedtuple


k_max = 7
peak_x = 0.17
peak_y = 0.3

problem = He5
bases = [mom,]
basis_size = 45

ks =[[],[]]
problem.V0 = -47.
#He5.Vso = -9.58

contour = triangle_contour(peak_x, peak_y, k_max, basis_size/3)

QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))
for m,basis in enumerate(bases):
    print basis.name
    eigvals, eigvecs = basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        ks[m].append(sp.sqrt(2*problem.mass*eigvals[i])    )
    
    plt.plot(sp.real(ks[m]),sp.imag(ks[m]),'*')
    
    # if basis == osc:
    #     wavef = basis.gen_wavefunction(eigvecs[:,0])
    # else:
    #     wavef = mom.gen_wavefunction(eigvecs[:,0], contour, Q)
    # 
    # print "Lowest energy:", eigvals[0]
    # r = sp.linspace(0, 10, 100)
    # plt.title(basis.name)
    # plt.plot(r, abs(wavef(r)*r))

plt.show()