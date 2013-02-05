from imports import *
from scipy import integrate
import unittest

def norm(f, start, stop, weight = lambda x: 1):
    def integrand(x):
        return sp.conjugate( f(x) ) * f(x) * weight(x);
    (N, _) = sp.integrate.quad(integrand, start, stop)
    return N


def normalize(f, start, stop, weight = lambda x: 1):
    N = norm(f, start, stop, weight = weight)
    return lambda x: f(x) / sp.sqrt(N)
    
#
# Tests
#
   
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
    