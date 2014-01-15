# This Python file uses the following encoding: utf-8
import scipy as sp
from scipy import linalg, integrate

def matrix_from_function(function, order, dtype=complex, 
                         hermitian=False, symmetric=False):
    """
    Given a function of two integer indexes, generates a square
    matrix of a certain order using that function. If the matrix 
    is known to be hermitian or symmetric, this can be entered to
    almost double the speed of the algorithm.
    
    This function would greatly benefit from being run in parallel.
    """
    matrix = sp.empty((order, order), dtype)
    for i in xrange(order):
        limit = (i + 1) if hermitian or symmetric else order
        for j in xrange(limit):
            matrix[i, j] = function(i, j)
            if hermitian:
                matrix[j, i] = sp.conj(matrix[i, j])
            elif symmetric:
                matrix[j, i] = matrix[i, j]
    return matrix
    
def symmetric(matrix):
    """
    Returns true if the matrix is symmetric.
    """
    return sp.allclose(matrix, matrix.T)
    
def hermitian(matrix):
    """
    Returns true if the matrix is hermitian.
    """
    return sp.allclose(matrix, sp.conj(matrix.T))

def eigensolve(H, hermitian=False):
    """
    Given a hamiltonian matrix, calculates energies and 
    eigenvectors (≈ wavefunctions) and sorts them by 
    the real part of the energy in ascending order.
    """
    if hermitian:
        energies, eigenvectors = linalg.eigh(H)
    else:
        energies, eigenvectors = linalg.eig(H)
        indexes = energies.argsort()
        energies = sp.real_if_close(energies[indexes])
        eigenvectors = eigenvectors[:, indexes]
    return energies, eigenvectors