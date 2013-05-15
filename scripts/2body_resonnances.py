from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq, energies
from nhqm.bases.gen_contour import gauss_contour, triangle_contour
from collections import namedtuple

k_max = 7

problem = He5
basis = mom
basis_size = 10
bs = basis_size
tail_size = 20
momenta =[[],[]]
eigvals = [[],[]]
eigvecs = [[],[]]
corner_x = 0.5
corner_y = -0.45j
problem.V0 = -47.
#He5.Vso = -9.58

p12res_exp = 0.35 -0.2j
p32res_exp = 0.17 -0.03j

def find_closest(vec, val):
    #returns index, value
    index= min(range(len(vec)), key=lambda i: abs(vec[i]-val))
    return index, vec[index]

def calc():

            
    contour = gauss_contour([0,corner_y, corner_x + corner_y, corner_x, k_max], [bs,bs,bs,tail_size] )
    QNums = namedtuple('qnums', 'l j J M E m')
    Qs = [QNums(l=1, j=0.5, J=0, M=0, 
              m=[-1.5, -0.5, 0.5, 1.5], 
              E=range(basis_size)),
        QNums(l=1, j=1.5, J=0, M=0, 
                  m=[-1.5, -0.5, 0.5, 1.5], 
                  E=range(basis_size))]
                    
    for m,Q in enumerate(Qs):
        eigvals[m], eigvecs[m] = basis.solution(contour, problem, Q)
        for i in xrange(len(eigvals[m])):
            momenta[m].append(sp.sqrt(2 * problem.mass * eigvals[m][i]) )
            if sp.imag(momenta[m][i]) > 0:
                momenta[m][i] = sp.conj(momenta[m][i])
    return eigvals, eigvecs, momenta    
        
def export_resonance():
    _, eigvecs, momenta = calc()
    i32, res32 = find_closest(momenta[1], p32res_exp)
    i12, res12 = find_closest(momenta[0], p12res_exp)
    return eigvecs[0][i12],eigvecs[1][i32]    
        
if __name__ == '__main__':        
    _, _, momenta = calc()
    i32, res32 = find_closest(momenta[1], p32res_exp)
    i12, res12 = find_closest(momenta[0], p12res_exp)
    print "p3/2 ressonance:",res32
    print "p1/2 resonnance:",res12

    plt.plot(sp.real(momenta[1]),sp.imag(momenta[1]),'b*')
    plt.plot(sp.real(momenta[1][i32]),sp.imag(momenta[1][i32]),'r*')
    plt.show()