from __future__ import division
import scipy as sp
from scipy import linalg, integrate
from scipy.special.specfun import csphjy
from scipy.misc import factorial
from nhqm.helpers import memoize

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

def j_l(l, x):
    """Spherical bessel."""
    L = 1 if l == 0 else l
    _, j_l, _, _, _ = csphjy(L, x)
    return j_l[l]

# Limit on number of iterations for computing C-G coefficients
CG_LIMIT = 50

@memoize
def clebsch_gordan(j1, j2, m1, m2, J, M):
    """Computes the Clebsch-Gordan coefficient
    <j1 j2; m1 m2|j1 j2; J M>.

    For reference see
    http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients."""
    if M != m1 + m2 or not abs(j1 - j2) <= J <= j1 + j2:
        return 0
    c1 = sp.sqrt((2*J+1) * factorial(J+j1-j2) * factorial(J-j1+j2) * \
              factorial(j1+j2-J)/factorial(j1+j2+J+1))
    c2 = sp.sqrt(factorial(J+M) * factorial(J-M) * factorial(j1-m1) * \
              factorial(j1+m1) * factorial(j2-m2) * factorial(j2+m2))
    c3 = 0.
    for k in range(CG_LIMIT):
        use = True
        d = [0, 0, 0, 0, 0]
        d[0] = j1 + j2 - J - k
        d[1] = j1 - m1 - k
        d[2] = j2 + m2 - k
        d[3] = J - j2 + m1 + k
        d[4] = J - j1 -m2 + k
        prod = factorial(k)
        for arg in d:
            if arg < 0:
                use = False
                break
            prod *= factorial(arg)
        if use:
            c3 += (-1)**k/prod
    return c1*c2*c3
    
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
#   Tests
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