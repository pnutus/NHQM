from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.QM_helpers import matrix_from_function, j_l, energies, berggren_norm

name = "MomSpace"

integration_order = 60
integration_range = 20

def solution(contour, problem, Q):
    H = hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    norm_eigvecs = sp.empty(sp.shape(eigvecs), complex)
    for col in xrange(len(eigvecs)):
         norm_eigvecs[:,col] = eigvecs[:,col] / berggren_norm(eigvecs[:,col])
    return eigvals, norm_eigvecs

def hamiltonian(contour, problem, Q):
    points, weights = contour
    def H_func(i, j):
        return H_element(points[i], points[j], weights[i], weights[j], problem, Q)
    return matrix_from_function(H_func, order=len(points), symmetric=True)
    
def H_element(k, k_prim, weight, weight_prim, problem, Q):
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    integral, _ = fixed_quad(integrand, 0, integration_range, n = integration_order,
                                args=(k, k_prim, V, Q.l, Q.j))
    return diagonal + 2*k*k_prim*sp.sqrt(weight*weight_prim) * integral / sp.pi

@sp.vectorize
def integrand(r, k, k_prim, V, l, j):
    return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)

def gen_wavefunction(eigvec, contour, Q):
    points, weights = contour
    @sp.vectorize
    def wavefunction(r):
        j_ls = [j_l(Q.l, points[n]*r) for n in xrange(len(points))]
        # A factor 1j**Q.l is removed for convenience
        return sp.sqrt(2/sp.pi)*sp.sum(sp.sqrt(weights)*points*j_ls*eigvec)
    return wavefunction
