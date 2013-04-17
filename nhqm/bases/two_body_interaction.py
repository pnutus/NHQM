from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom

V0 = -1e9 #MeV
beta = 1 # fm^-2

def interaction(a, b, c, d):
    return V0*(sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c])

def gen_matrix(eigvecs, contour, Q, verbose=False):
    order = len(eigvecs)
    points, weights = contour
    
    def V_func(i, j):
        return V_sep(points[i], points[j], Q)
    V_matrix = matrix_from_function(V_func, order)
    
    global sep_M
    def sep_M_func(i, j):
        return V_matrix.dot(eigvecs[i]*weights*points**2) \
                       .dot(eigvecs[j]*weights*points**2)
    sep_M = 2 / sp.pi * matrix_from_function(sep_M_func, order)
    return sep_M
    
def matrix_from_function(function, order, dtype=complex):
    matrix = sp.empty((order, order), dtype)
    for i in xrange(order):
        for j in xrange(order):
            matrix[i, j] = function(i, j)
    return matrix
    
def V_sep(k, k_prim, Q):
    args = (k, k_prim, potential, Q.l, Q.j)
    integral, _ = fixed_quad(mom.integrand, 0, 10, n = 20, args=args)
    return integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
