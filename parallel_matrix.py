from imports import *
from IMap import IMap
import multiprocessing as mp
import unittest

def one(i,j):
    return i*j + 1

def above_diagonal(row, col):
    return (col - row) > 0
    
def matrix_indices(k, matrix_dim):
    """Generates matrix indices from a single number."""
    return (k // matrix_dim, k % matrix_dim) 
    
def mat_range(matrix_dim, existing_dim, hermitian = False):
    """Generates indices for """
    
    for i in xrange(matrix_dim**2):
        if ((( i % matrix_dim ) < existing_dim )  and (( i // matrix_dim ) < existing_dim )):
            yield None
        elif hermitian and above_diagonal(*matrix_indices(i, matrix_dim)):
            yield None
        else:
            yield i
        
def hermitize(M):
    """Mirrors a upper/lower triangular matrix into an hermitian."""
    return M - sp.diag(M.A.diagonal()) + M.H

def element_from_index(k, f, matrix_dim):
    if k is None:
        return 0
    else:
        return f(*matrix_indices(k, matrix_dim))
    
def parallel_matrix(f, matrix_dim, hermitian = False, existing_dim = 0 ):
    """Calculates a matrix in parallel using multiple processes."""
    
            
    # p = mp.Pool()
    #    result = p.imap(element_from_index, mat_range(matrix_dim))
    imap = IMap()
    result = imap.add_task(element_from_index, \
                mat_range(matrix_dim, existing_dim, hermitian), (f, matrix_dim))
    imap.start()
    
    M = sp.mat( sp.empty((matrix_dim, matrix_dim)) )
    for k, val in enumerate(result):
        i, j = matrix_indices(k, matrix_dim)
        M[i, j] = val
    if hermitian:
        M = hermitize(M)
    return M

#
# Tests
#

liten = 3
stor = 7
mat = parallel_matrix (one, stor, False, liten)
print mat
   
"""class hermitizeTests(unittest.TestCase):
    
    def setUp(self):
        self.M = sp.mat("[1, 2, 3; 4, 5, 6; 7, 8, 9]")
    
    def testLower(self):
        L = sp.tril(self.M)
        H = hermitize(L)
        A = sp.mat("[1, 4, 7; 4, 5, 8; 7, 8, 9]")
        self.assertTrue( sp.array_equal(H, A) )
        
    def testUpper(self):
        U = sp.triu(self.M)
        H = hermitize(U)
        A = sp.mat("[1, 2, 3; 2, 5, 6; 3, 6, 9]")
        self.assertTrue( sp.array_equal(H, A) )
        
    def testComplex(self):
        pass
        
class parallel_matrixTests(unittest.TestCase):
    
    def setUp(self):
        self.correct = sp.mat("0, 1, 2; 1, 2, 3; 2, 3, 4")
        
    def testParallel(self):
        M = parallel_matrix(simple, 3)
        self.assertTrue( sp.array_equal(M, self.correct) )
        
    def testHermitian(self):
        M = parallel_matrix(simple, 3, hermitian = True)
        self.assertTrue( sp.array_equal(M, self.correct) )
        
def simple(i, j):
    return i + j
        
if __name__ == '__main__':
    unittest.main()"""
    