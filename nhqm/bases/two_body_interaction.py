from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function

V0 = -200 #MeV
beta = 1 # fm^-2

def gen_matrix(eigvecs, Q, basis, contour=None):
    order = len(eigvecs)
    
    if basis == mom:
        points, weights = contour
        dk = sp.sqrt(weights) * points
        def V_func(i, j):
            integral, _ = fixed_quad(mom.integrand, 0, 10, n = 60,
                        args=(points[i], points[j], potential, Q.l, Q.j))
            return 2 / sp.pi * integral
    elif basis == osc:
        dk = 1
        def V_func(i, j):
            integral, _ = fixed_quad(osc.integrand, 0, 10, n = 30,
                                        args=(i, j, potential, Q))
            return integral
    V_ij = matrix_from_function(V_func, order, symmetric=True)
    
    def sep_M_func(i, j):
        return (eigvecs[:,i]*dk).dot(V_ij).dot(eigvecs[:,j]*dk)
        
    sep_M = matrix_from_function(sep_M_func, order, symmetric=True)
    return sp.sqrt(V0) * sep_M
    
@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
