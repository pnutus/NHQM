from imports import *
from scipy import integrate, optimize, linalg
from sympy.physics.sho import R_nl
import sympy as sym
import math

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

@sp.vectorize
def V(r):
    return - 1 / r
    
def Kdelta(x, y):
    return x == y

@sp.vectorize
def H_element(n, n_prim, l, s, eps):
    """Returns matrix element of the Hamiltonian"""
    #n' => m; n => n
    # http://docs.sympy.org/0.7.2/modules/physics/sho.html
    rsq = eps**2 * .5 * ( \
        sp.sqrt(n * (n + l + .5)) * Kdelta(n_prim, n - 1) + \
        2.*(2*n + l + 1.5) * Kdelta(n_prim, n) + \
        sp.sqrt( (n + 1) * (n + l + 1.5) ) * Kdelta(n_prim, n + 1) \
        )
    def integrand(r):
        return R_nl(n_prim, l, 1, r) * V(r) * R_nl(n, l, 1, r) * r**2.
    
    (integral, _) = integrate.quad(integrand, 0, sp.inf)
    return rsq + sp.sqrt(2) * eps * integral

@sp.vectorize
def ground_energy(omega):
    return H_element(0, 0, 0, .5, omega)

# Plots the ground energy as a function of the frequency parameter
# x = sp.linspace(0,2)
# plt.plot(x, ground_energy(x))
# plt.show()

# Varies the frequency to minimize the ground energy
# res = optimize.minimize_scalar(ground_energy, bounds=(0, 2), method='bounded')
# print res.x # 0.752252778064

minomega = 0.752252778064

def hamiltonian(n):
    """Calculates the n x n hamiltonian matrix."""
    n = 10
    H = sp.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            H[i, j] = H[j, i] = H_element(i, j, 0, .5, minomega)
    return H

# Diagonalize
H = hamiltonian(3)
energies = linalg.eigvals(H)
print H
print energies

# Save to file
sp.savetxt("matris.txt",H)

# Profit???