from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function

name = "Gaussian"



V0 = 100 # MeV
r0 = 2 # fm

def gen_interaction(sp_states):
    global dk
    points, weights = sp_states[0].contour
    dk = sp.sqrt(weights) * points
    sep_M = gen_matrix(sp_states)
    def interaction(a, b, c, d):
        return sep_M[a.id, c.id] * sep_M[b.id, d.id] + sep_M[a.id, d.id] * sep_M[b.id, c.id]
    return interaction

def gen_matrix(sp_states):
    order = len(sp_states[0].eigvec)
    def V_func(i, j):
        if sp_states[i].basis == mom:
            points, weights = sp_states[i].contour
            dk = sp.sqrt(weights) * points
            integral, _ = fixed_quad(mom.integrand, 0, 20, n = 60,
                        args=(points[i], points[j], potential, sp_states[i].l, sp_states[i].j))
            return 2 / sp.pi * integral
        elif sp_states[i].basis == osc:
            #dk = 1
            integral, _ = fixed_quad(osc.integrand, 0, 20, n = 60,
                                        args=(i, j, potential, Q))
            return integral
    V_ij = matrix_from_function(V_func, order, symmetric=True)
    
    def sep_M_func(i, j):
        def phi(index):
            return sp_states[index].eigvec * dk
        return phi(i).dot(V_ij).dot(phi(j))    
        
    def sep_M_func_sum(i, j):
        """this function is just for comparison, doesn't work with HO, can't be arsed to fix it """
        res = []
        points, weights = sp_states[i].contour
        for k in xrange(len(sp_states[i].eigvec)):
            for l in xrange(len(sp_states[i].eigvec)):
                temp = sp.sqrt(weights[k]) * points[k] * sp_states[i].eigvec[k]
                temp *= sp.sqrt(weights[l]) * points[l] * sp_states[j].eigvec[l]
                res.append( temp * V_ij[k,l])
        return sp.sum(res)        
        
    sep_M = matrix_from_function(sep_M_func, order, symmetric=True)
    return sp.sqrt(-V0) * sep_M
    
@sp.vectorize
def potential(r, l, j):
    return sp.exp(-(r/r0)**2)
    
