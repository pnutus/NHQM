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
    integral, _ = fixed_quad(integrand, 0, 10, n = 20,
                                args=(k, k_prim, V, Q.l, Q.j))
    return diagonal + 2 * k_prim**2 * weight / sp.pi * integral

@sp.vectorize
def integrand(r, k, k_prim, V, l, j):
    return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)

def gen_basis_function(problem, l = 0, j = .5):
    # TODO
    # contour/complex
    # Not normalized!
    # Doesn't work with G-L quadrature
    def basis_function(r, k, weight):
        return 4*sp.pi/(2*sp.pi)**1.5 * k**2 * weight*j_l(l, k*r)
    return basis_function