from imports import *
import HO_basis

@sp.vectorize
def V(r, l, s):
    return - 1. / r

H = HO_basis.hamiltonian(V, 10, verbose=True)
E = HO_basis.energies(H)
print "H =", H
print "E =", E
print "GT = ", E[0] * 2.