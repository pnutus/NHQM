from imports import *
from scipy import integrate
import math
import unittest

def norm(f, start, stop, weight = lambda x: 1):
    def integrand(x):
        return sp.conjugate( f(x) ) * f(x) * weight(x);
    (N, _) = sp.integrate.quad(integrand, start, stop)
    return N


def normalize(f, start, stop, weight = lambda x: 1):
    N = norm(f, start, stop, weight = weight)
    return lambda x: f(x) / math.sqrt(N)
    
#
# Tests
#
   
class SinTests(unittest.TestCase):
    
    def setUp(self):
        self.a = 0
        self.b = 5
        self.f = lambda x: math.sin(math.pi / self.b * x)
    
    def testNorm(self):
        N = norm(self.f, self.a, self.b)
        self.failUnless(N == 2.5)
    
    def testNormalize(self):
        g = normalize(self.f, self.a, self.b)
        N = norm(g, self.a, self.b)
        self.failUnless(N == 1.0)

if __name__ == '__main__':
    unittest.main()
    