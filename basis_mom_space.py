from imports import *
from scipy import integrate, linalg, special
import config
from scipy.special import sph_jn

@sp.vectorize
def fourier_integrand(r, p_0, V, l, j):
    if p_0 == 0:
        return V(r, l, j) * r**2
    else:
        return r / p_0 * V(r, l, j) * sp.sin(p_0 * r)

@sp.vectorize
def cos_integrand(x, p, p_prim, V, l, j):
    p_0 = sp.sqrt(p**2 + p_prim**2 - 2*p*p_prim*x)
    V_tilde, _ = integrate.fixed_quad(fourier_integrand, 0, 5, n = 20, \
                 args = (p_0, V, l, j))
    return V_tilde * special.legendre(l)(x)

@sp.vectorize
def l_zero_integrand(r, p, p_prim, V):
    return V(r, 0, .5) * sp.sin(p * r) * sp.sin(p_prim * r)

def H_element(n, n_prim, problem, l = 0, j = .5, step_size = 0.25):
    p, p_prim = [(x + 1)*step_size for x in (n, n_prim)]
    diagonal = p**2 / (2 * problem.mass) * (n == n_prim)
    V = problem.potential
    if l == 0: # faster, only one integral
        integral, _ = integrate.fixed_quad(l_zero_integrand, 0, 10, n = 20, \
                        args = (p, p_prim, V))
        integral *=  2 / (p * p_prim)
    else:
        integral, _ = integrate.fixed_quad(cos_integrand, -1, 1, n = l + 5, \
                args = (p, p_prim, V, l, j))
    return diagonal + p_prim**2 * step_size / sp.pi * integral

def gen_basis_function(problem, step_size = .25, l = 0, j = .5):
    def basis_function(r, n):
        p_n = (n + 1) * step_size
        jn, _ = sph_jn(l, p_n*r)
        return 4*sp.pi / (2*sp.pi)**1.5 * p_n**2 * step_size * jn[-1]
    return basis_function

