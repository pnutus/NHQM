from imports import *
import HO_basis

@sp.vectorize
def V(r, l, s):
    return - 1. / r
    
# @sp.vectorize
# def ground_energy(eps):
#     """First matrix element as a function of epsilon"""
#     return HO_basis.H_element2(0, 0, 0, .5, eps, V)
# 
# x =  sp.linspace(0, 4)
# plt.plot(x, ground_energy(x))
# plt.show()

H = HO_basis.hamiltonian(V, 10, verbose=True)
E = HO_basis.energies(H)
print "H =", H
print "E =", E
print "GT = ", E[0] * 2.