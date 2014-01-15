# This Python file uses the following encoding: utf-8
from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.helpers.matrix import matrix_from_function
from nhqm.helpers.quantum import j_l
from scipy.special.orthogonal import p_roots

class MomentumBasis:
    name = "Momentum Basis"
    hermitian = False
    
    def __init__(self, contour, int_order=60, int_range=20):
        self.contour = contour
        self.integration_order = int_order # Number of points in integration
        self.integration_range = int_range #fm
    
    def hamiltonian(self, problem, quantum_numbers):
        points, weights = self.contour
        def H_func(i, j):
            return self.H_element(
                points[i], points[j], weights[i], weights[j], 
                problem.mass, problem.potential, quantum_numbers
            )
        return matrix_from_function(H_func, order=len(points), symmetric=True)
    
    def H_element(self, k, k_prim, weight, weight_prim, mass, V, Q):
        diagonal = k**2 / (2 * mass) * (k == k_prim)
        integral, _ = fixed_quad(integrand, 0, self.integration_range, 
                                    n = self.integration_order,
                                    args=(k, k_prim, V, Q.l, Q.j))
        return diagonal + 2*k*k_prim*sp.sqrt(weight*weight_prim) * integral / sp.pi
    
    def normalize_eigenvectors(self, eigvecs):
        normed_eigvecs = sp.empty(sp.shape(eigvecs), complex)
        for col in xrange(len(eigvecs)):
             normed_eigvecs[:,col] = berggren_normalize(eigvecs[:, col])
        return normed_eigvecs

    def gen_wavefunction(self, eigvec, Q):
        points, weights = self.contour
        @sp.vectorize
        def wavefunction(r):
            j_ls = [j_l(Q.l, points[n]*r) for n in xrange(len(points))]
            # A complex phase factor 1j**Q.l is removed for convenience
            # as it has no physical significance
            return sp.sqrt(2/sp.pi)*sp.sum(sp.sqrt(weights)*points*j_ls*eigvec)
        return wavefunction

@sp.vectorize
def integrand(r, k, k_prim, V, l, j):
    return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)

#
### Contour
#
# A contour is a collection of points and corresponding weights. 
#

def gauss_contour(vertices, points_per_segment):
    """
    On the line connecting each pair of vertices, generates
    points and weights using Gauss-Legendre quadrature.
    """
    segment_count = len(vertices) - 1
    points = weights = sp.empty(0, complex)
    for i in range(segment_count):
        try:    (x, w) = p_roots(points_per_segment[i])
        except: (x, w) = p_roots(points_per_segment)
        a = vertices[i]
        b = vertices[i + 1]
        scaled_x = (x * (b - a) + (a + b))/2
        scaled_w = w * (b - a)/2
        points = sp.hstack((points, scaled_x))
        weights = sp.hstack((weights, scaled_w))
    return (points, weights)
        
def triangle_contour(peak, k_max, points_per_segment):
    """
    Generates a triangular contour with a specified peak 
    ending at k_max, with a certain number of points per segment.
    
          _ _ _ _ k_max
    \    /
     \  /
      \/
     peak
    
    """
    vertices = [0, peak, 2 * sp.real(peak), k_max]
    return gauss_contour(vertices, points_per_segment)

#
### Berggren
#
# Peculiarities of the Berggren basis
#
    
def berggren_norm(x):
    """
    Computes the norm of a vector using the non-conjugated
    scalar product. 
    < x | y > = ∑ x_i y_i
    < x | y > ≠ ∑ x_i y_i^*
    """
    return sp.sqrt(sp.dot(x, x))
