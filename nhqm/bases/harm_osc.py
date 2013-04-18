from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.optimize import minimize_scalar
from scipy.misc import factorial, factorial2
from scipy.special import genlaguerre
from collections import namedtuple
from nhqm.helpers import matrix_from_function

name = "HarmOsc"

QNums = namedtuple('qnums', 'l j n')


def hamiltonian(order, problem, Q):
    eps = 1.06384608107
    def H_func(i, j):
        return H_element_mosh(i, j, problem, eps, Q)
    return matrix_from_function(H_func, order, hermitian=True)

def H_element_mosh(n, n_prim, problem, eps, Q):
    l = Q.l
    rsq = eps**2 / 2 * (
        sp.sqrt(n * (n + l + .5)) * (n_prim == (n - 1))
        + (2*n + l + 1.5) * (n_prim == n)
        + sp.sqrt( (n + 1) * (n + l + 1.5) ) * (n_prim == (n + 1))
        )
    integrand = gen_integrand(n, n_prim, .5, problem.potential, Q)
    (integral, _) = fixed_quad(integrand, 0, 10, n = 20)
    return rsq + sp.sqrt(2) * eps * integral

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
    integrand = gen_integrand(n, n_prim, nu, V, Q)
    (integral, _) = fixed_quad(integrand, 0, 10, n = 20)
    
    return omega / 2 * result + float(integral)

def gen_integrand(n, n_prim, nu, V, Q):
    R_nl = gen_R_nl(n, Q.l, nu)
    R_nl_prim = gen_R_nl(n_prim, Q.l, nu)
    def integrand(r):
        return R_nl(r)*V(r, Q.l, Q.j)*R_nl_prim(r)*r**2
    return integrand

def gen_R_nl(n, l, nu):
    laguerre = genlaguerre(n, l + .5)
    def R_nl(r):
        norm = sp.sqrt(
                        sp.sqrt(2*nu**3/sp.pi)*
                        2**(n + 2*l + 3)*factorial(n, exact=True)*
                        nu**l/factorial2(2*n + 2*l + 1, exact=True)
                      )
        lagge = laguerre(2*nu*r**2)
        return norm * r**l * sp.e**(-nu * r**2) * lagge
    return R_nl


def gen_basis_function(problem, omega, l = 0, j = .5):
    nu = problem.mass * omega / 2
    def basis_function(r, n):
        return R_nl(n, l, nu, r)
    return basis_function