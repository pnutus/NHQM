from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function


V0 = -200 # MeV
r0 = 2 # fm

def gen_interaction(eigvecs, Q, sp_basis, contour=None):
    sep_M = gen_matrix(eigvecs, Q, sp_basis, contour)
    def interaction(a, b, c, d):
        return sep_M[a, c] * sep_M[b, d]
    return interaction

def gen_matrix(eigvecs, Q, sp_basis, contour=None):
    order = len(eigvecs)
    
    if sp_basis == mom:
        points, weights = contour
        dk = sp.sqrt(weights) * points
        def V_func(i, j):
            integral, _ = fixed_quad(mom.integrand, 0, 20, n = 60,
                        args=(points[i], points[j], potential, Q.l, Q.j))
            return 2 / sp.pi * integral
    elif sp_basis == osc:
        dk = 1
        def V_func(i, j):
            integral, _ = fixed_quad(osc.integrand, 0, 20, n = 60,
                                        args=(i, j, potential, Q))
            return integral
    V_ij = matrix_from_function(V_func, order, symmetric=True)
    
    
    def sep_M_func(i, j):
        def phi(index):
            return eigvecs[:,index] * dk
        return phi(i).dot( V_ij.dot(phi(j)) )
        
        
    def sep_M_func_sum(i, j):
        """this function is just for comparison, doesn't work with HO, can't be arsed to fix it """
        res = []
        points, weights = contour
        for k in xrange(len(eigvecs[:,i])):
            for l in xrange(len(eigvecs[:,j])):
                temp = sp.sqrt(weights[k]) * points[k] * eigvecs[k,i]
                temp *= sp.sqrt(weights[l]) * points[l] * eigvecs[l,j]
                res.append( temp * V_ij[k,l])
        return sp.sum(res)        
        
    sep_M = matrix_from_function(sep_M_func, order, symmetric=True)
    return sp.sqrt(V0) * sep_M
    
@sp.vectorize
def potential(r, l, j):
    return sp.exp(-(r/r0)**2)
    
