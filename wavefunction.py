from imports import *

def gen_wavefunction(eigvec, basis_function):
    length = len(eigvec)
    @sp.vectorize
    def wavefunction(r):
        result = sp.empty(length)
        for (n, weight) in enumerate(eigvec):
            result[n] = weight * basis_function(r, n)
        return sp.sum(result)
    return wavefunction