from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.special.specfun import csphjy
import numpy as np

name = "MomSpace"

def j_l(l, k):
    """Spherical bessel."""
    if l == 0:
    l1 = 1 if l == 0 else l
    else:
        l1 = l
    _, j_l, _, _, _ = csphjy(l1, k) 
    return j_l[l]

@sp.vectorize
def integrand(r, k, k_prim, V, l, j):
    return r**2 * j_l(l, sp.conj(k)*r) \
    return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)
    
def H_element_contour(k, k_prim, weight, problem, l = 0, j = .5):
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    integral, _ = fixed_quad(integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j))
    return diagonal + 2 * k_prim**2 * weight / sp.pi * integral

def H_element(n, n_prim, step, *args, **kwargs):
    k, k_prim = [(x + 1)*step for x in (n, n_prim)]
    return H_element_contour(k, k_prim, step, *args, **kwargs)

def gen_basis_function(problem, l = 0, j = .5):
    # TODO
    # contour/complex
    # Not normalized!
    # Doesn't work with G-L quadrature
    def basis_function(r, k, weight):
        return 4*sp.pi/(2*sp.pi)**1.5 * k**2 * weight*j_l(l, k*r)
    return basis_function


# Berggren contour:    

@sp.vectorize
def triangle(x, peak_x, peak_y, k_max):
    if x < 2*peak_x:
        return peak_y * (abs(x/peak_x - 1) - 1)
    else:
        return 0

def gen_simple_contour(peak_x, peak_y, k_max, points):
    xs = sp.linspace(0, k_max, points + 1)[1:]
    ys = triangle(xs, peak_x, peak_y, k_max)
    return xs + 1j * ys
