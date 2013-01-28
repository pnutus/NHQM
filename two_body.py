from imports import *
from math import *
import scipy.integrate
from sympy.physics.sho import R_nl

def Kdelta(x, y):
    return x == y



def Helement(n, m, l, s):
    #n' => m; n => n
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    def integrand(r):
        return np.conjugate(R_nl(m,l,1,r)) * \
            V(r,l,s) * R_nl(n,l,1,r) * r**2
    Hho = (2*n + l + 1.5) * Kdelta(m, n)
    integral = sp.integrate.quad(integrand, 1e-2, 1)
    return rsq + Hho + integral

print Helement(0, 0, 0, .5)