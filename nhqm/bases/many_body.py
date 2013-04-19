from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, combinations_with_replacement
from collections import namedtuple
from nhqm.bases.fermion_state import FermionState
import nhqm.bases.two_body_interaction as n_n
from nhqm.QM_helpers import clebsch_gordan
            


# THIS IS FOR FERMIONS

QNums = namedtuple('qnums', 'l J M j m k')
SP = namedtuple('sp', ['m', 'k'])

def hamiltonian(quantum_numbers, 
                eigvals, eigvecs, contour, 
                num_particles=2, verbose=False):
    
    if verbose: print "Generating matrix of separable interactions"
    n_n.gen_matrix(eigvecs, contour, quantum_numbers)
    
    # tuples of k-indexes: (3, 6), (4, 4)
    mb_k = list(combinations_with_replacement(quantum_numbers.k, 
                                              num_particles))
    order = len(mb_k)
    
    if verbose: print "Generating many-body H"
    H = sp.zeros((order, order), complex)
    for i, k_bra in enumerate(mb_k):
        for j, k_ket in enumerate(mb_k):
            H[i,j] = coupled_H_elem(k_bra, k_ket, eigvals, quantum_numbers)
    return H

def coupled_H_elem(k_bra, k_ket, eigvals, Q):
    k1, k2 = k_bra
    k1_prim, k2_prim = k_ket
    
    result = 0
    for m in Q.m:
        bra = FermionState([SP(m, k1), SP(Q.M - m, k2)])
        if bra.sign == 0: continue
        
        for m_prim in Q.m:
            ket = FermionState([SP(m_prim, k1_prim), 
                                SP(Q.M - m_prim, k2_prim)])
            if ket.sign == 0: continue
            temp = clebsch_gordan(Q.j, Q.j, m, Q.M - m, Q.J, Q.M)
            temp *= clebsch_gordan(Q.j, Q.j, m_prim, Q.M - m_prim, Q.J, Q.M)
            temp *= H_elem(bra, ket, eigvals)
            result += temp
    return result
            
def H_elem(bra, ket, eigvals):
    # H_1, one-body interaction
    one_body = 0
    if bra == ket:
        for sp_state in bra:
            one_body += eigvals[sp_state.k]
    
    # H_2, two-body interaction
    two_body = sum(sign * n_n.interaction(a, b, c, d)
                    for (a, b, c, d, sign) in two_body_indexes(bra, ket))       
    return one_body + two_body



def two_body_indexes(bra, ket):
    result = []
    if len(set(bra) - set(ket)) > 2:
        return []
    for annihilated in combinations(ket, 2):
        for created in combinations(bra, 2):
            new_ket = (ket
                        .annihilate(annihilated)
                        .create(created))
            if new_ket.states == bra.states:
                sign = new_ket.sign
                a,b = created
                c,d = annihilated
                result.append( (a,b,c,d,sign) )
    return result
    
import unittest

class RedTests(unittest.TestCase):
    
    def setUp(self):
        self.bra32 = FermionState([1,2,3]) 
        self.ket32 = FermionState([1,3,4])
        self.res32 = [(1, 2, 1, 4, -1), (2, 3, 3, 4, 1)]
        
        self.bra31 = FermionState([1,2,3]) 
        self.ket31 = FermionState([3,5,6])
        self.res31 = [(1,2, 5,6, +1)]

        self.bra33 = FermionState([1,2,3]) 
        self.ket33 = FermionState([1,2,3])
        self.res33 = [(1,2, 1,2, +1), (1,3,1,3,1), (2,3,2,3,1)]
        
        
    def test32(self):
        res = two_body_indexes(self.bra32, self.ket32)
        
        self.assertEquals(res, self.res32 )
        
    def test31(self):
        res = two_body_indexes(self.bra31, self.ket31)
        
        self.assertEquals(res, self.res31 ) 
        
    def test33(self):
        res = two_body_indexes(self.bra33, self.ket33)
        
        self.assertEquals(res, self.res33 )       
        

            
if __name__ == '__main__':
    print "kaptenkvant 4 lyfe"
    unittest.main()       

