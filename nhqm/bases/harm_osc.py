from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.misc import factorial, factorial2
from scipy.special import genlaguerre
from collections import namedtuple
from nhqm.QM_helpers import matrix_from_function, energies

name = "HarmOsc"

QNums = namedtuple('qnums', 'l j n')

def solution(order, problem, Q):
    H = hamiltonian(order, problem, Q)
    eigvals, eigvecs = energies(H)
    
    # Don't repeat this?
    omega = problem.HO_omega
    nu = problem.mass * omega / 2
    R_nls = [gen_R_nl(n, Q.l, nu) for n in xrange(order)]
    
    wavefunctions = [gen_wavefunction(eigvec, R_nls) for eigvec in eigvecs]
    return eigvals, wavefunctions

def gen_wavefunction(eigvec, R_nls):
    def wavefunction(r):
        return sum(R_nls[n](r)*eigvec[n] for n in xrange(len(eigvec)))
    return wavefunction

def hamiltonian(order, problem, Q):
    omega = problem.HO_omega
    nu = problem.mass * omega / 2
    R_nls = [gen_R_nl(n, Q.l, nu) for n in xrange(order)]
    def H_func(i, j):
        return H_element(i, j, problem, omega, Q, R_nls)
    return matrix_from_function(H_func, order, hermitian=True)
    
def H_element(n, n_prim, problem, omega, Q, R_nls):
    if n == n_prim:
        result = (2*n + Q.l + 1.5)
    elif n_prim == (n - 1):
        result = sp.sqrt( n*(n + Q.l + .5) )
    elif n_prim == (n + 1):
        result = sp.sqrt( (n + 1)*(n + Q.l + 1.5) )
    else:
        result = 0

    V = problem.potential
    (integral, _) = fixed_quad(integrand, 0, 10, n = 20, 
                                args=(n, n_prim, V, Q, R_nls))
    return omega / 2 * result + float(integral)
    
def integrand(r, n, n_prim, V, Q, R_nls):
    return R_nls[n](r)*V(r, Q.l, Q.j)*R_nls[n_prim](r)*r**2
    
def gen_R_nl(n, l, nu):
    norm = sp.sqrt(
                    sp.sqrt(2*nu**3/sp.pi)*
                    2**(n + 2*l + 3)*factorial(n, exact=True)*
                    nu**l/factorial2(2*n + 2*l + 1, exact=True)
                  )
    laguerre = genlaguerre(n, l + .5)
    def R_nl(r):
        return norm * r**l * sp.exp(-nu * r**2) * laguerre(2*nu*r**2)
    return R_nl


def gen_basis_function(problem, Q):
    omega = problem.HO_omega
    nu = problem.mass * omega / 2
    R_nls = [gen_R_nl(n, Q.l, nu) for n in xrange(order)]
    def basis_function(r, n):
        return R_nls[n](r)
    return basis_function
