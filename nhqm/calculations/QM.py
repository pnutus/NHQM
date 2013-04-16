from __future__ import division
import scipy as sp
from scipy import linalg, integrate
from scipy.special.orthogonal import p_roots

def hamiltonian(element_function, args=None, 
                order=20, hermitian=False):
    """Given an element_function(n, n_prim) calculates energies."""
    H = sp.empty((order, order), complex)
    for n in xrange(order):
        limit = (n + 1) if hermitian else order
        for n_prim in xrange(limit):
            H[n, n_prim] = element_function(n, n_prim, *args)
            if hermitian:
                H[n_prim, n] = H[n, n_prim]
    return H

def contour_hamiltonian(element_function, contour, args=None):
    zip_contour = zip(*contour)
    order = len(zip_contour)
    H = sp.empty((order, order), complex)
    for (n, (k, _)) in enumerate(zip_contour):
        for (n_prim, (k_prim, w)) in enumerate(zip_contour):
            H[n, n_prim] = element_function(k, k_prim, w, *args)
    return H



def energies(H, hermitian=False):
    """
    Given a hamiltonian matrix, calculates energies and 
    sorts them by their real part.
    """
    if hermitian:
        eigvals, eigvecs = linalg.eigh(H)
    else:
        eigvals, eigvecs = linalg.eig(H)
        indexes = eigvals.argsort()
        eigvals = sp.real_if_close(eigvals[indexes])
        eigvecs = eigvecs[:, indexes]
    return eigvals, eigvecs

def gen_wavefunction(eigvec, basis_function, contour=None):
    length = len(eigvec)
    if contour == None:
        @sp.vectorize
        def wavefunction(r):
            iterable = (basis_function(r, n) 
                        for n in range(length))
            basis_vec = sp.fromiter(iterable, complex, count=length)
            return sp.dot(basis_vec, eigvec)
    else:
        @sp.vectorize
        def wavefunction(r):
            iterable = (basis_function(r, k, w)
                            for (k, w) in zip(*contour))
            basis_vec = sp.fromiter(iterable, complex, count=length)
            return sp.dot(basis_vec, eigvec)
    return wavefunction
    

# Berggren contour:    

@sp.vectorize
def triangle(x, peak_x, peak_y, k_max):
    if x < 2*peak_x:
        return peak_y * (abs(x/peak_x - 1) - 1)
    else:
        return 0

def naive_triangle_contour(peak_x, peak_y, k_max, points):
    xs = sp.linspace(0, k_max, points + 1)[1:]
    ys = triangle(xs, peak_x, peak_y, k_max)
    points = xs + 1j * ys
    weights = calculate_steps(points)
    return (points, weights)
    
def calculate_steps(contour):
    shifted = sp.roll(contour, 1)
    shifted[0] = 0
    return contour - shifted
    
def gauss_contour(vertices, order):
    """
    Generates a contour along the line segments between
    vertices using Gauss-Legendre quadrature.
    """
    (x, w) = p_roots(order)
    num_segments = len(vertices) - 1
    points = weights = sp.empty(0, complex)
    for i in range(num_segments):
        a = vertices[i]
        b = vertices[i + 1]
        scaled_x = (x * (b - a) + (a + b))/2
        scaled_w = w * (b - a)/2
        points = sp.hstack((points, scaled_x))
        weights = sp.hstack((weights, scaled_w))
    return (points, weights)
    
def triangle_contour(peak_x, peak_y, k_max, order):
    vertices = [0, peak_x - 1j*peak_y, 2*peak_x, k_max]
    return gauss_contour(vertices, order)

def absq(x):
    return sp.real(x)**2 + sp.imag(x)**2

def norm(f, start, stop, weight = lambda x: 1):
    def integrand(x):
        return absq(f(x)) * weight(x);
    (N, _) = sp.integrate.quad(integrand, start, stop)
    return N


def normalize(f, start, stop, weight = lambda x: 1):
    N = norm(f, start, stop, weight = weight)
    return lambda x: f(x) / sp.sqrt(N)
    
#
# Tests
#
   
import unittest
   
class SinTests(unittest.TestCase):
    
    def setUp(self):
        self.a = 0
        self.b = 5
        self.f = lambda x: sp.sin(sp.pi / self.b * x)
    
    def testNorm(self):
        N = norm(self.f, self.a, self.b)
        self.assertEquals(N, (self.b - self.a) * .5 )
    
    def testNormalize(self):
        g = normalize(self.f, self.a, self.b)
        N = norm(g, self.a, self.b)
        self.assertEquals(N, 1.0)
        
class ComplexTests(unittest.TestCase):
    pass
    
class WeightTests(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()
    