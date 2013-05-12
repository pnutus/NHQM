from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import harm_osc as osc
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

def gen_matrix(eigvecs, Q, verbose=False):
    order = len(eigvecs)
    
    def V_func(i, j):
        return V_sep(i, j, Q)
    V_matrix = matrix_from_function(V_func, order, symmetric=True)
    
    global sep_M
    def sep_M_func(i, j):
        dk = 1
        return V_matrix.dot(eigvecs[j]*dk).dot(eigvecs[i]*dk)
        
    sep_M = 2 / sp.pi * matrix_from_function(sep_M_func, order, 
                                             symmetric=True)
    return sep_M
    
def V_sep(n, n_prim, Q):
    V = potential
    integral, _ = fixed_quad(osc.integrand, 0, 10, n = 60,
                                args=(n, n_prim, V, Q))
    return integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
