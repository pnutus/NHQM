from __future__ import division
import scipy as sp
from scipy import linalg, integrate


def matrix_from_function(function, order, dtype=complex, 
                         hermitian=False, symmetric=False):
    matrix = sp.empty((order, order), dtype)
    for i in xrange(order):
        limit = (i + 1) if hermitian or symmetric else order
        for j in xrange(limit):
            matrix[i, j] = function(i, j)
            if hermitian:
                matrix[j, i] = sp.conj(matrix[i, j])
            elif symmetric:
                matrix[j, i] = matrix[i, j]
    return matrix
    
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
    
#
#   Decorators
#
    
import collections
import functools

class memoize(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)