from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.helpers import matrix_from_function
from nhqm.QM_helpers import j_l
from collections import namedtuple

name = "MomSpace"

QNums = namedtuple('qnums', 'l j k')

def hamiltonian(contour, problem, Q):
    points, weights = contour
    def H_func(i, j):
        return H_element(points[i], points[j], weights[j], problem, Q)
    return matrix_from_function(H_func, order=len(points))
    
def H_element(k, k_prim, weight, problem, Q):
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    
    V = problem.potential
    def integrand(r):
        return r**2 * j_l(Q.l, k*r) * j_l(Q.l, k_prim*r) * V(r, Q.l, Q.j)
    integral, _ = fixed_quad(integrand, 0, 10, n = 20)
    
    return diagonal + 2 * k_prim**2 * weight / sp.pi * integral

def gen_basis_function(problem, l = 0, j = .5):
    # TODO
    # contour/complex
    # Not normalized!
    # Doesn't work with G-L quadrature
    def basis_function(r, k, weight):
        return 4*sp.pi/(2*sp.pi)**1.5 * k**2 * weight*j_l(l, k*r)
    return basis_function