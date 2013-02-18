from imports import *
from scipy import linalg

def hamiltonian(element_function, args=None, order=20, hermitian=False):
    """Given an element_function(n, n_prim) calculates energies."""
    H = sp.empty((order, order))
    for n in xrange(order):
        limit = (n + 1) if hermitian else order
        for n_prim in xrange(limit):
            H[n, n_prim] = element_function(n, n_prim, *args)
            if hermitian:
                H[n_prim, n] = H[n, n_prim]
    return H

def energies(H, hermitian=False):
    """Given a hamiltonian matrix, calculates energies and sorts them."""
    if hermitian:
        eigvals, eigvecs = linalg.eigh(H)
    else:
        eigvals, eigvecs = linalg.eig(H)
        indexes = sp.real_if_close(eigvals).argsort()
        eigvals = sp.real_if_close(eigvals[indexes])
        eigvecs = eigvecs[:, indexes]
    return eigvals, eigvecs

