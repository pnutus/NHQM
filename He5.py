from imports import *
import HO_basis
import config

config.mass = 0.019272
config.omega_interval = (0, 50)

@sp.vectorize
def V(r, l, s):
    """Woods-Saxon potential for the He5 nucleus"""
    V0 = -70 # MeV
    Vso = -7.5 # MeV
    r0 = 2. # fm
    d = .65 # fm
    f = 1/(1 + sp.e**((r - r0)/d))
    j = abs(l + s)
    spin_orbit = .5 * ( j*(j + 1) - l*(l + 1) - .75 )
    return f * (V0 - 4 * Vso * spin_orbit * (f - 1) / (d * r))

def plotV():
    x = sp.linspace(0.1, 10, 100)
    plt.plot(x, V(x, 0, 0.5), 'k', x, V(x, 1, .5), 'r', x, V(x, 1, -.5), 'b')
    # plt.axis([7, 10, -2, 2])
    plt.show()
    
def calc_energies():
    H = HO_basis.hamiltonian(V, 10, l = 1, s = .5, verbose=True)
    E = HO_basis.energies(H)
    # print "H =", H
    print "E =", E


if __name__ == '__main__':
    pass
    calc_energies()
    # plotV()
    
    

