from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.special import sph_jn
import numpy as np

name = "MomSpace"
    
@sp.vectorize
def integrand(r, k, k_prim, V, l, j):
    j_l, _ = sph_jn(l, sp.conj(k)*r) 
    j_l_prim, _ = sph_jn(l, k_prim*r)
    return sp.real(r**2 * j_l[-1] * j_l_prim[-1] * V(r, l, j))
    
def H_element_contour(k, k_prim, step, problem, l = 0, j = .5):
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    integral, _ = fixed_quad(integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j))
    return diagonal + 2 * k_prim**2 * step / sp.pi * integral

def H_element(n, n_prim, step, *args, **kwargs):
    k, k_prim = [(x + 1)*step for x in (n, n_prim)]
    return H_element_contour(k, k_prim, step, *args, **kwargs)

def gen_basis_function(step, problem, l = 0, j = .5):
    # TODO
    # contour/complex
    # Not normalized!
    def basis_function(r, n):
        k_n = (n + 1) * step
        jn, _ = sph_jn(l, k_n*r)
        return 4*sp.pi / (2*sp.pi)**1.5 * k_n**2 * step * jn[-1]
    return basis_function