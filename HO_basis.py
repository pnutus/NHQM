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
    rsq = eps**2 * .5 * ( \
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
def H_element(n, n_prim, l, j, omega, V):
    """Returns matrix element of the Hamiltonian"""
    result = 0
    if n == n_prim:
        result += (2*n + l + 1.5)
    elif n_prim == (n - 1):
        result += sp.sqrt( n*(n + l + .5) )
    elif n_prim == (n + 1):
        result += sp.sqrt( (n + 1)*(n + l + 1.5) )
    
    # hbar = 1
    # mass = 1
    nu = config.mass * omega * .5
    @sp.vectorize
    def integrand(r):
        return R_nl(n_prim, l, nu, r) * V(r, l, j) * R_nl(n, l, nu, r) * r**2.
    (integral, _) = integrate.fixed_quad(integrand, 0, 10, n = 20)
    # (integral, _) = integrate.quad(integrand, 0, sp.inf)
    return omega / 2. * result + float(integral)


def optimized_eps(V, l, j):
    """Returns the epsilon that gives the minimum ground energy"""
    @sp.vectorize
    def ground_energy(eps):
        """First matrix element as a function of epsilon"""
        return H_element(0, 0, l, j, eps, V)
    res = optimize.minimize_scalar(ground_energy, \
                method = 'bounded', bounds = config.omega_interval)
    return float(res.x) # 1.06384608107

def hamiltonian(V, N, l = 0, j = .5, verbose=False):
    """Calculates the N x N hamiltonian matrix."""
    eps = optimized_eps(V, l, j)
    if verbose:
        print "omega", eps
    H = sp.empty((N, N))
    for i in range(N):
        if verbose:
            print "rad", i + 1
        for k in range(i+1):
            H[i, k] = H[k, i] = H_element(i, k, l, j, eps, V)
    return H
    
def energies(H):
    return linalg.eigh(H, eigvals_only = True)
    
def ground_state_wavefunction(H, V, l, j):
    _, eigvecs = linalg.eigh(H)
    eigvec = eigvecs[:, 0]
    dim = H.shape[0]
    omega = optimized_eps(V, l, j)
    nu = config.mass * omega * .5
    @sp.vectorize
    def wavefunction(r):
        funcs = sp.empty(dim)
        for i in xrange(dim):
            funcs[i] = R_nl(i, l, nu, r)
        return sp.sum(eigvec * funcs)
    return wavefunction
    
if __name__ == '__main__':
    # import H_moshinsky
    import He5