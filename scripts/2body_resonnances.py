from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import absq, energies
from nhqm.bases.contours import gauss_contour, triangle_contour
from collections import namedtuple


k_max = 7


problem = He5
basis = mom
basis_size = 15
bs = basis_size
tail_size = 25
ks =[[],[]]
problem.V0 = -47
#He5.Vso = -9.58

corner_x = 0.5
corner_y = -0.45j
problem.V0 = -47.
#He5.Vso = -9.58

p12res_exp = 0.35 -0.2j
p32res_exp = 0.17 -0.03j

def find_closest(vec, val):
    #returns index, value
    index = min(range(len(vec)), key=lambda i: abs(vec[i]-val))
    return index, vec[index]
            
contour = gauss_contour([0,corner_y, corner_x + corner_y, corner_x, k_max], [bs,bs,bs,tail_size] )
QNums = namedtuple('qnums', 'l j J M E m')
Qs = [QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size)),
    QNums(l=1, j=0.5, J=0, M=0, 
              m=[-1.5, -0.5, 0.5, 1.5], 
              E=range(basis_size))]
                    
for m,Q in enumerate(Qs):
    print "l:", Q.l, "j:", Q.j
    
    eigvals, eigvecs = basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        ks[m].append(sp.sqrt(2*problem.mass*eigvals[i]))
        if sp.imag(ks[m][i]) > 0:
            ks[m][i] = sp.conj(ks[m][i])
        
i32, res32 = find_closest(ks[0], p32res_exp)
i12, res12 = find_closest(ks[1], p12res_exp)
print "p3/2 ressonance:",res32
print "p1/2 resonnance:",res12

plt.plot(sp.real(ks[1]),sp.imag(ks[1]),'b*')
plt.plot(sp.real(ks[1][i12]),sp.imag(ks[1][i12]),'r*')
#plt.show()