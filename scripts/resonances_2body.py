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
corner_x = 0.5
corner_y = -0.45j
problem.V0 = -47.
default_contour = gauss_contour([0,corner_y, corner_x + corner_y, corner_x, k_max], [bs,bs,bs,tail_size] )
#He5.Vso = -9.58

#redundant quantum numbers?
QNums = namedtuple('qnums', 'l j J M E m')
Q12 = QNums(l=1, j=0.5, J=0, M=0, 
    m=[-1.5, -0.5, 0.5, 1.5], 
    E=range(basis_size))
    
Q32 = QNums(l=1, j=1.5, J=0, M=0, 
    m=[-1.5, -0.5, 0.5, 1.5], 
    E=range(basis_size))    

p12res_exp = 0.35 -0.2j
p32res_exp = 0.17 -0.03j

def find_closest(vec, val):
    #returns index, value
    index= min(range(len(vec)), key=lambda i: abs(vec[i]-val))
    
    """
    TODO:
    Raise a flag or w/e if the closest value is too close
    i e no pole was found
    """
    return index, vec[index]

def calc(Q, contour = default_contour):
    momenta =[]
    eigvals, eigvecs = basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        momenta.append(sp.sqrt(2 * problem.mass * eigvals[i]) )
        if sp.imag(momenta[i]) > 0:
            momenta[i] = sp.conj(momenta[i])
            
    return eigvals, eigvecs, momenta    
        
def export_resonance_p12(contour=default_contour):
    eigvals, eigvecs, momenta = calc(Q12, contour)
    i12, res12 = find_closest(momenta, p12res_exp)
    return eigvals[i12], eigvecs[i12]
        
if __name__ == '__main__':        
    # _, _, momenta12 = calc(Q12)
    # _, _, momenta32 = calc(Q32) 
    # i32, res32 = find_closest(momenta32, p32res_exp)
    # i12, res12 = find_closest(momenta12, p12res_exp)
    # print "p3/2 ressonance:",res32
    # print "p1/2 resonnance:",res12
    print "export p1/2:", export_resonance_p12()
    # 
    # plt.plot(sp.real(momenta32),sp.imag(momenta32),'b*')
    # plt.plot(sp.real(momenta32[i32]),sp.imag(momenta32[i32]),'r*')
    # plt.show()