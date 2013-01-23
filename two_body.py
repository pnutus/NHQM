from imports import *
from math import *
import sympy
import scipy.integrate
from sympy.physics.sho import R_nl

@np.vectorize
def V(r, l, s):
    V0 = -70e6
    Vso = -7.5e6
    r0 = 2e-15
    d = .65e-15
    g = e**((r - r0)/d)
    f = 1/(1 + g)
    df = - g/( d * (1 + g)**2 )
    return V0*f - 4 * Vso * l * s * df / r

def Kdelta(x, y):
    if x==y:
        return 1
    else:
        return 0

x = np.linspace(1e-16, 5e-15, 1000)
plt.plot(x, V(x, 1, .5))
plt.axis([0, 5e-15, -1e30, 0])
plt.show()

def Helement(n, nn, l):
    #n' => nn; n => n
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    def integrand(r, n, nn, l):
        return np.conjugate(R_nl(nn,l,1,r))*V(r)*R_nl(n,l,1,r)*r**2
    rsq = sqrt(n*(n+l+0.5))*Kdelta(nn, n-1) + sqrt((n+1)*(n+l+1.5))*Kdelta(nn,n+1)
    Hho = (2*n+l+1.5)*Kdelta(nn,n)
    integral = integrate.quad(integrand, start, end)
    return rsq + Hho + integral

