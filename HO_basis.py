from imports import *
from scipy import integrate, optimize, linalg
from sympy.physics.sho import R_nl

def Kdelta(x, y):
    return x == y

@sp.vectorize
def H_element_mosh(n, n_prim, l, s, eps, V):
    """Returns matrix element of the Hamiltonian"""
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    rsq = eps**2 * .5 * ( \
        sp.sqrt(n * (n + l + .5)) * Kdelta(n_prim, n - 1) + \
        (2*n + l + 1.5) * Kdelta(n_prim, n) + \
        sp.sqrt( (n + 1) * (n + l + 1.5) ) * Kdelta(n_prim, n + 1) \
        )
    @sp.vectorize
    def integrand(r):
        return R_nl(n_prim, l, .5, r) * V(r, l, s) * R_nl(n, l, .5, r) * r**2.
    
    (integral, _) = integrate.fixed_quad(integrand, 0, 10, n = 20)
    return rsq + sp.sqrt(2) * eps * integral

@sp.vectorize
def H_element(n, n_prim, l, s, omega, V):
    """Returns matrix element of the Hamiltonian"""
    result = 0
    if n == n_prim:
        result += (2*n + l + 1.5)
    elif n_prim == (n - 1):
        result += sp.sqrt( n*(n + l + .5) )
    elif n_prim == (n + 1):
        result += sp.sqrt( (n + 1)*(n + l + 1.5) )
    
    # hbar = 1
    mass = 1 
    nu = mass * omega * .5
    @sp.vectorize
    def integrand(r):
        return R_nl(n_prim, l, nu, r) * V(r, l, s) * R_nl(n, l, nu, r) * r**2.
    (integral, _) = integrate.fixed_quad(integrand, 0, 10, n = 20)
    return omega / 2. * result + float(integral)


def optimized_eps(V, l, s):
    """Returns the epsilon that gives the minimum ground energy"""
    @sp.vectorize
    def ground_energy(eps):
        """First matrix element as a function of epsilon"""
        return H_element(0, 0, l, s, eps, V)
    res = optimize.minimize_scalar(ground_energy, bounds=(0, 2), method='bounded')
    return float(res.x) # 1.06384608107

def hamiltonian(V, N, l = 0, s = .5, verbose=False):
    """Calculates the N x N hamiltonian matrix."""
    eps = optimized_eps(V, l, s)
    if verbose:
        print "omega", eps
    H = sp.empty((N, N))
    for i in range(N):
        if verbose:
            print "rad", i+1
        for j in range(i+1):
            H[i, j] = H[j, i] = H_element(i, j, l, s, eps, V)
    return H


def energies(H):
    return linalg.eigh(H, eigvals_only = True)
    
def ground_state_wavefunction(H):
    """docstring for ground_state_wavefunction"""
    pass

if __name__ == '__main__':
    import H_moshinsky