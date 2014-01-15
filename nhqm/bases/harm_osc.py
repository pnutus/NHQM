from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.misc import factorial, factorial2
from scipy.special import genlaguerre
from nhqm.helpers.matrix import matrix_from_function
from nhqm.helpers.decorators import memoize

class HarmOscBasis:
    
    name = "Harmonic Oscillator Basis"
    hermitian = True
    
    def __init__(self, 
                 basis_state_count, mass, omega,
                 int_order=70, int_range=20):
        self.basis_state_count = basis_state_count 
        self.mass = mass
        self.omega = omega
        self.integration_order = int_order # Number of points in integration
        self.integration_range = int_range
        
    
    def hamiltonian(self, problem, Q):
        """
        Given a basis size and a sph.symm. problem,
        generates a Hamiltonian matrix in the Spherical
        Harmonic Oscillator basis.
        """
        nu = self.mass * self.omega / 2
        R_nls = [gen_R_nl(n, Q.l, nu) for n in xrange(self.basis_state_count)]
        def H_func(i, j):
            return self.H_element(i, j, problem.potential, 
                                  problem.HO_omega, Q, R_nls)
        return matrix_from_function(H_func, self.basis_state_count, 
                                    hermitian=True)
    
    def H_element(self, n, n_prim, V, omega, Q, R_nls):
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
        (integral, _) = fixed_quad(integrand, 0, 
                                    self.integration_range, 
                                    n = self.integration_order, 
                                    args=(n, n_prim, V, Q, R_nls))
        return omega / 2 * result + integral

    def gen_wavefunction(self, eigvec, Q):
        """
        Given a solution eigenvector, generates plottable
        wavefunctions psi(r).
        """
        ns = range(len(eigvec))
        nu = self.mass * self.omega / 2
        R_nls = [gen_R_nl(n, Q.l, nu) for n in ns]
        def wavefunction(r):
            return sum( R_nls[n](r)*eigvec[n] for n in ns )
        return wavefunction

    
    
    
    
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
