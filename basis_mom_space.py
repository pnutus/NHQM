from imports import *
from scipy import integrate, linalg, special
from scipy.special import sph_jn

name = "MomSpace"
default_step_size = 0.25
    
@sp.vectorize
def integrand(r, p, p_prim, V, l, j):
    j_l, _ = sph_jn(l, p*r)
    j_l_prim, _ = sph_jn(l, p_prim*r)
    return r**2 * j_l[-1] * j_l_prim[-1] * V(r, l, j)
    

def H_element(n, n_prim, problem, l = 0, j = .5, step_size = default_step_size):
    p, p_prim = [(x + 1)*step_size for x in (n, n_prim)]
    diagonal = p**2 / (2 * problem.mass) * (n == n_prim)
    V = problem.potential
    integral, _ = integrate.fixed_quad(integrand, 0, 10, \
                                n = 20, args=(p, p_prim, V, l, j))
    return diagonal + 2 * p_prim**2 * step_size / sp.pi * integral

def gen_basis_function(problem, step_size = default_step_size,\
                            l = 0, j = .5):
    # Felnormaliserad!
    def basis_function(r, n):
        p_n = (n + 1) * step_size
        jn, _ = sph_jn(l, p_n*r)
        return 4*sp.pi / (2*sp.pi)**1.5 * p_n**2 * step_size * jn[-1]
    return basis_function

