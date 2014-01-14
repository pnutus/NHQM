from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.misc import factorial, factorial2
from scipy.special import genlaguerre
from nhqm.QM_helpers import matrix_from_function, energies
from nhqm.helpers import memoize

name = "HarmOsc"
integration_order = 70 # Number of points used in integration
integration_range = 20 # fm

def solution(basis_size, problem, Q):
    """
    Given a basis size and a sph.symm. problem,
    returns the eigenvalues and corresponding
    eigenvectors for the Hamiltonian in the 
    Spherical Harmonic Oscillator basis.
    """
    try: 
        basis_size = len(basis_size[0])
    except: 
        pass
    H = hamiltonian(basis_size, problem, Q)
    eigvals, eigvecs = energies(H, hermitian=True)
    return eigvals, eigvecs

def gen_wavefunction(eigvec, Q, problem, contour=None):
    """
    Given a solution eigenvector, generates plottable
    wavefunctions psi(r).
    """
    ns = range(len(eigvec))
    nu = problem.mass * problem.HO_omega / 2
    R_nls = [gen_R_nl(n, Q.l, nu) for n in ns]
    def wavefunction(r):
        return sum( R_nls[n](r)*eigvec[n] for n in ns )
    return wavefunction

def hamiltonian(basis_size, problem, Q):
    """
    Given a basis size and a sph.symm. problem,
    generates a Hamiltonian matrix in the Spherical
    Harmonic Oscillator basis.
    """
    nu = problem.mass * problem.HO_omega / 2
    R_nls = [gen_R_nl(n, Q.l, nu) for n in xrange(basis_size)]
    def H_func(i, j):
        return H_element(i, j, problem, problem.HO_omega, Q, R_nls)
    return matrix_from_function(H_func, basis_size, hermitian=True)
    
def H_element(n, n_prim, problem, omega, Q, R_nls):
    """
    Calculates the value of one element in the hamiltonian matrix.
    """
    if n == n_prim:
        result = (2*n + Q.l + 1.5)
    elif n_prim == (n - 1):
        result = sp.sqrt( n*(n + Q.l + .5) )
    elif n_prim == (n + 1):
        result = sp.sqrt( (n + 1)*(n + Q.l + 1.5) )
    else:
        result = 0
    
    def integrand(r, n, n_prim, V, Q, R_nls):
        return sp.conj(R_nls[n](r))*V(r, Q.l, Q.j)*R_nls[n_prim](r)*r**2
    
    V = problem.potential
    (integral, _) = fixed_quad(integrand, 0, 
                                integration_range, 
                                n = integration_order, 
                                args=(n, n_prim, V, Q, R_nls))
    return omega / 2 * result + integral
    
@memoize
def gen_R_nl(n, l, nu):
    """
    Generates a radial eigenfunction for the Spherical
    Harmonic Oscillator. Used to speed up calculations.
    """
    norm = sp.sqrt(
                    sp.sqrt(2*nu**3/sp.pi)*
                    2**(n + 2*l + 3)*factorial(n, exact=True)*
                    nu**l/factorial2(2*n + 2*l + 1, exact=True)
                  )
    laguerre = genlaguerre(n, l + .5)
    def R_nl(r):
        return norm * r**l * sp.exp(-nu * r**2) * laguerre(2*nu*r**2)
    return R_nl
