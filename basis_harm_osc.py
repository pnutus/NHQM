from __future__ import division
from imports import *
from scipy import integrate, optimize, linalg
from sympy.physics.sho import R_nl
import config

def Kdelta(x, y):
    return x == y

@sp.vectorize
def H_element_mosh(n, n_prim, l, j, eps, V):
    """Returns matrix element of the Hamiltonian"""
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    rsq = eps**2 / 2 * ( \
        sp.sqrt(n * (n + l + .5)) * Kdelta(n_prim, n - 1) + \
        (2*n + l + 1.5) * Kdelta(n_prim, n) + \
        sp.sqrt( (n + 1) * (n + l + 1.5) ) * Kdelta(n_prim, n + 1) \
        )
    @sp.vectorize
    def integrand(r):
        return R_nl(n_prim, l, .5, r) * V(r, l, j) * R_nl(n, l, .5, r) * r**2.
    
    (integral, _) = integrate.fixed_quad(integrand, 0, 10, n = 20)
    return rsq + sp.sqrt(2) * eps * integral

@sp.vectorize
def H_element(n, n_prim, problem, omega, l = 0, j = .5):
    """Returns matrix element of the Hamiltonian"""
    
    if n == n_prim:
        result = (2*n + l + 1.5)
    elif n_prim == (n - 1):
        result = sp.sqrt( n*(n + l + .5) )
    elif n_prim == (n + 1):
        result = sp.sqrt( (n + 1)*(n + l + 1.5) )
    else:
        result = 0

    nu = problem.mass * omega / 2
    V = problem.potential
    @sp.vectorize
    def integrand(r):
        return R_nl(n_prim, l, nu, r) * V(r, l, j) * R_nl(n, l, nu, r) * r**2.
    (integral, _) = integrate.fixed_quad(integrand, 0, 10, n = 20)
    
    return omega / 2 * result + float(integral)


def optimal_osc_freq(problem, l = 0, j = .5):
    """Returns the epsilon that gives the minimum ground energy"""
    @sp.vectorize
    def ground_energy(omega):
        return H_element(0, 0, problem, omega, l, j)
    res = optimize.minimize_scalar(ground_energy, \
                method = 'bounded', bounds = problem.omega_interval)
    return float(res.x) # 1.06384608107

def gen_basis_function(problem, omega, l = 0, j = .5):
    nu = problem.mass * omega / 2
    def basis_function(r, n):
        return R_nl(n, l, nu, r)
    return basis_function