from imports import *
from central_problem import CentralProblem

@sp.vectorize
def V(r, l, j):
    """Woods-Saxon potential for the He5 nucleus"""
    V0 = -70 # MeV
    Vso = -7.5 # MeV
    r0 = 2. # fm
    d = .65 # fm
    f = 1/(1 + sp.e**((r - r0)/d))
    spin_orbit = .5*(j*(j + 1) - l*(l + 1) - .75)
    return f * (V0 - 4*Vso*spin_orbit*(f - 1) / (d * r))
    
problem = CentralProblem()
problem.potential = V
problem.mass = 0.019272
problem.eV_factor = 1e6
problem.omega_interval = (0, 50)

from scipy import fftpack
import HO_basis

def plotV():
    x = sp.linspace(0.1, 10, 100)
    plt.plot(x, V(x, 0, 0.5), 'k', x, V(x, 1, 1.5), 'r', x, V(x, 1, .5), 'b')
    # plt.axis([7, 10, -2, 2])
    plt.show()
    
def calc_energies():
    l = 0
    j = .5
    H = HO_basis.hamiltonian(V, 10, l, j, verbose=True)
    E = HO_basis.energies(H)
    # print "H =", H
    print "E =", E
    wavef = HO_basis.ground_state_wavefunction(H, V, l, j)
    r = sp.linspace(0, 10, 1024)
    y = fftpack.fft(wavef(r))
    plt.plot(y[1:128])
    plt.show()
    
    # plt.plot(r, wavef(r) )
    # plt.show()


if __name__ == '__main__':
    pass
    # calc_energies()
    # plotV()
    
    

