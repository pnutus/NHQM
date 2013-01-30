from imports import *
import HO_basis

@sp.vectorize
def V(r, l, s):
    """Woods-Saxon potential for the He5 nucleus"""
    V0 = -70 # MeV
    Vso = -7.5 # MeV
    r0 = 2. # fm
    d = .65 # fm
    g = e**((r - r0)/d)
    f = 1/(1 + g)
    df = 
    j = l + s
    spin_orbit = .5 * ( j*(j + 1) - l*(l + 1) - s*(s + 1) )
    return f * (V0 - 4 * Vso * spin_orbit / (d * r))

# x = sp.linspace(1e-16, 5e-15, 100)
# plt.plot(x, V(x, 1, .5), 'r', x, V(x, 0, .5), 'k')
# plt.axis([0, 5e-15, -1e37, 0])
# plt.figure(2)
# plt.plot(x, V(x, 0, .5), 'k')
# plt.show()