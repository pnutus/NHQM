from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import He5, H_atom
from nhqm.QM_helpers import absq, energies
from nhqm.bases.gen_contour import gauss_contour, triangle_contour
from collections import namedtuple
from scipy.optimize import fmin
import numpy
import matplotlib.pyplot as plt

k_max = 7

problem = He5
basis = mom
basis_size = 15*3
bs = basis_size
problem.V0 = -47.
peak_x = 0.35
peak_y = 0.4
default_contour = triangle_contour(peak_x, peak_y, k_max, basis_size/3)


QNums = namedtuple('qnums', 'l j J M E m')
Q12 = QNums(l=1, j=0.5, J=0, M=0, 
    m=[-1.5, -0.5, 0.5, 1.5], 
    E=range(basis_size))
    
Q32 = QNums(l=1, j=1.5, J=0, M=0, 
    m=[-1.5, -0.5, 0.5, 1.5], 
    E=range(basis_size))    

#energi
Ep12res_exp = 4.9-0.5j
Ep31res_exp = 0.89-0.3j
#momentum
#p12res_exp = 0.35 -0.2j
#p32res_exp = 0.17 -0.03j
p12res_exp = .435-.022
p32res_exp = .1878-.031

def find_closest(vec, val):
    #returns index, value
    index= min(range(len(vec)), key=lambda i: abs(vec[i]-val))
    return index, vec[index]

def solve_2b(contour,basis,problem,Q):
    if basis == osc:
        H = osc.hamiltonian(basis_size, problem, Q)
    else:
        H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    return eigvals, eigvecs

def calc(Q, contour = default_contour):
    momenta =[]
    eigvals, eigvecs = solve_2b(contour,basis, problem, Q)
    for i in xrange(len(eigvals)):
        momenta.append(sp.sqrt(2 * problem.mass * eigvals[i]) )
        if sp.imag(momenta[i]) > 0:
            momenta[i] = sp.conj(momenta[i])
            
    return eigvals, eigvecs, momenta

def export_resonance_p12(contour=default_contour):
    eigvals, eigvecs, momenta = calc(Q12, contour)
    i12, _ = find_closest(momenta, p12res_exp)
    return momenta[i12]

def export_resonance_p32(contour=default_contour):
    eigvals, eigvecs, momenta = calc(Q32, contour)
    i32, _ = find_closest(momenta, p32res_exp)
    return momenta[i32]

resonance=p12res_exp

def f_to_min(x):
    problem.V0=x[0]
    problem.Vso=x[1]
    momenta=export_resonance_p12()
    # ret = ev-resonance
    
    ret=numpy.absolute(momenta-resonance)
    print ret
    return ret

def minimize():
    return fmin(f_to_min,[-47.,-7.5],xtol=0.0001, ftol=0.0001, maxiter=100)
# print minimize()
# ev,_=export_resonance_p32()
problem.V0=-51.30070425
problem.Vso=-12.36063252
ev=export_resonance_p12()
print ev

# ev,_,_=calc(Q32, default_contour)
# ev=sp.sqrt(2*problem.mass*ev)
# plt.plot(sp.real(ev),sp.imag(ev),'*')
# plt.show()
