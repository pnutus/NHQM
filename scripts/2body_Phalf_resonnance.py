from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq, energies
from nhqm.bases.gen_contour import gauss_contour, triangle_contour
from collections import namedtuple

problem = He5
bases = [mom,]
basis_size = 45
k_max = 7
ks =[[],[]]
peak_x = 0.3
peak_y = -0.45j
problem.V0 = -47.
He5.Vso = -9.58

contour = gauss_contour([0,peak_x +peak_y, 2*peak_x, k_max], [30,30,40] )

QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))
                    
          
for m, basis in enumerate(bases):
    print basis.name
    eigvals, eigvecs = basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        ks[m].append(sp.sqrt(2*problem.mass*eigvals[i])    )
        #print i, ks[m][i]
        #42 (0.350565321605-0.200277314032j)
    plt.plot(sp.real(ks[m]),sp.imag(ks[m]),'*')
    res = res_index(eigvals)


    
    # if basis == osc:
#         wavef = basis.gen_wavefunction(eigvecs[:,0])
#     else:
#         wavef = mom.gen_wavefunction(eigvecs[:,0], contour, Q)
# 
#     print "Lowest energy:", eigvals[0]
#     r = sp.linspace(0, 10, 100)
#     plt.title(basis.name)
#     plt.plot(r, abs(wavef(r)*r))

plt.show()