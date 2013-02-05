from imports import *
from scipy import integrate, optimize, linalg
import config

config.mass = 0.019272

def index_to_coord(i, size):
    if not 0 <= i < size**3:
        raise IndexError('Index is out of bounds')
    x = i % size
    tmp = i // size
    y = tmp % size
    z = tmp / size
    return sp.array([x, y, z]) - (size - 1)/2.
    
def fourier_integrand(r, V, p, p_prim):
    p_0 = sp.linalg.norm(p - p_prim)
    if p_0 == 0:
        return 1 / (2 * sp.pi**2) * V(r, 0, 0) * r**2
    else:
        return 2 / (p_0 * (2*sp.pi)**2) * r * V(r, 0, 0) * sp.sin(p_0 * r)

def H_element(n, n_prim, size, V, step_size = .5):
    p = index_to_coord(n, size) * step_size
    p_prim = index_to_coord(n_prim, size) * step_size
    diagonal = sp.sum(p**2) / (2 * config.mass) * (n == n_prim)
    integral, _ = integrate.fixed_quad(fourier_integrand, 0, 10, n = 20, \
                 args = (V, p, p_prim))
    return diagonal + integral * step_size**3
    
# import He5

@sp.vectorize
def V(r, l, j):
    """Woods-Saxon potential for the He5 nucleus"""
    V0 = -70 # MeV
    Vso = -7.5 # MeV
    r0 = 2. # fm
    d = .65 # fm
    f = 1/(1 + sp.e**((r - r0)/d))
    spin_orbit = .5 * ( j*(j + 1) - l*(l + 1) - .75 )
    return f * (V0 - 4 * Vso * spin_orbit * (f - 1) / (d * r))

size = 10
H = sp.empty((size**3, size**3))
for i in xrange(size**3):
    print "rad", i
    for j in xrange(i + 1):
        H[i, j] = H[j, i] = H_element(i, j, size, V)

print sp.linalg.eigh(H, eigvals_only=True)[:3]


#
# Tests
#

import unittest
   
class IndexToCoordinateTests(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def testSizeOne(self):
        size = 1
        test = index_to_coord(0, size)
        correct = sp.array([0, 0, 0])
        self.assertTrue( sp.array_equal(test, correct) )
        
    def testSizeThree(self):
        size = 3
        test = index_to_coord(13, size)
        correct = sp.array([0, 0, 0])
        self.assertTrue( sp.array_equal(test, correct) )
        
    def testOverflow(self):
        with self.assertRaises(IndexError):
            size = 2
            test = index_to_coord(10, size)
        with self.assertRaises(IndexError):
            size = 5
            test = index_to_coord(-1, size)

# if __name__ == '__main__':
#     unittest.main()