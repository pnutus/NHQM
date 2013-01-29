
"""
Steps:
1. Find wavefunction of harm osc ground state as function of 
    frequency.
2. Calculate matrix element H_00 of ground state using that
    wavefunction.
3. Minimize the energy with respect to frequency.
4. Use that frequency with harm osc wavefunctions for excited
    states.
5. Use these to calculate more matrix elements of H.
"""


@np.vectorize
def V(r, l, s):
    """Woods-Saxon potential for the He5 nucleus"""
    V0 = -70e6
    Vso = -7.5e6
    r0 = 2e-15
    d = .65e-15
    g = e**((r - r0)/d)
    f = 1/(1 + g)
    df = - g/( d * (1 + g)**2 )
    return V0*f - 4 * Vso * l * s * df / r


# x = np.linspace(1e-16, 5e-15, 100)
# plt.plot(x, V(x, 1, .5), 'r', x, V(x, 0, .5), 'k')
# plt.axis([0, 5e-15, -1e37, 0])
# plt.figure(2)
# plt.plot(x, V(x, 0, .5), 'k')
# plt.show()

def Kdelta(x, y):
    return x == y

def integrand(r):
    return np.conjugate(R_nl(m,l,1,r)) * \
        V(r,l,s) * R_nl(n,l,1,r) * r**2

def Helement(n, m, l, s):
    #n' => m; n => n
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    Hho = (2*n + l + 1.5) * Kdelta(m, n)
    integral = sp.integrate.quad(integrand, 1e-2, 1)
    return Hho + integral

print Helement(0, 0, 0, .5)