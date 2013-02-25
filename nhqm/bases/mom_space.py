from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.special import sph_jn
import numpy as np

name = "MomSpace"
    
@sp.vectorize
def integrand(r, k, k_prim, V, l, j): #takes complex p, p_prim and returns real Re of integrand
    
    j_l, _ = sph_jn(l, k*r) 
    j_l_prim, _ = sph_jn(l, k_prim*r)
    return r**2 * j_l[-1] * j_l_prim[-1] * V(r, l, j)
    
def H_element(n, n_prim, problem, step_size, l = 0, j = .5):
    k, k_prim = [(x + 1)*step_size for x in (n, n_prim)] # discretized momentum space basis
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    integral, _ = fixed_quad(integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j))

    return diagonal + 2 * k_prim**2 * step_size / sp.pi * integral


def gen_basis_function(problem, step_size,\
                            l = 0, j = .5):
    # Felnormaliserad!
    def basis_function(r, n):
        k_n = (n + 1) * step_size
        jn, _ = sph_jn(l, k_n*r)
        return 4*sp.pi / (2*sp.pi)**1.5 * k_n**2 * step_size * jn[-1]
    return basis_function