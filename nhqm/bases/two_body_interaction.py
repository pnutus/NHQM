from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom
from nhqm.helpers import matrix_from_function
from nhqm.QM_helpers import j_l

V0 = -1e9 #MeV
beta = 1 # fm^-2

def interaction(a, b, c, d):
    return V0*(sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c])

def gen_matrix(eigvecs, contour, Q, verbose=False):
    order = len(eigvecs)
    points, weights = contour
    
    def V_func(i, j):
        return V_sep(points[i], points[j], Q)
    V_matrix = matrix_from_function(V_func, order, symmetric=True)
    
    global sep_M
    def sep_M_func(i, j):
        dk = weights * points**2
        return V_matrix.dot(eigvecs[i]*dk).dot(eigvecs[j]*dk)
    sep_M = 2 / sp.pi * matrix_from_function(sep_M_func, order, 
                                             symmetric=True)
    return sep_M
    
def V_sep(k, k_prim, Q):
    V = potential
    @sp.vectorize
    def integrand(r):
        return r**2 * j_l(Q.l, k*r) * j_l(Q.l, k_prim*r) * V(r, Q.l, Q.j)
    integral, _ = fixed_quad(integrand, 0, 10, n = 20)
    return integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
