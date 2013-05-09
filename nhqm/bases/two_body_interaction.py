from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom
from nhqm.QM_helpers import matrix_from_function, j_l

V0 = -200 #MeV
beta = 1 # fm^-2

def interaction(a, b, c, d):
    result = 0
    if a.m == c.m and b.m == d.m:
        result += sep_M[a.E, c.E]*sep_M[b.E, d.E]
    if a.m == d.m and b.m == c.m:
        result += - sep_M[a.E, d.E]*sep_M[b.E, c.E]
    return result*V0

def gen_matrix(eigvecs, contour, Q, verbose=False):
    order = len(eigvecs)
    points, weights = contour
    
    def V_func(i, j):
        return V_sep(points[i], points[j], Q)
    V_matrix = matrix_from_function(V_func, order, symmetric=True)
    
    global sep_M
    def sep_M_func(i, j):
        dk = sp.sqrt(weights) * points
        return V_matrix.dot(eigvecs[i]*dk).dot(eigvecs[j]*dk)
    sep_M = 2 / sp.pi * matrix_from_function(sep_M_func, order, 
                                             symmetric=True)
    return sep_M
    
def V_sep(k, k_prim, Q):
    V = potential
    integral, _ = fixed_quad(mom.integrand, 0, 10, n = 20,
                                args=(k, k_prim, V, Q.l, Q.j))
    return integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
