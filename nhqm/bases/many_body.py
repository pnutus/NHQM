from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, combinations_with_replacement
from collections import namedtuple
from nhqm.bases.fermion_state import FermionState
import nhqm.bases.two_body_interaction as n_n
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function
            


# THIS IS FOR FERMIONS

QNums = namedtuple('qnums', 'l J M j m E')
SP = namedtuple('sp', ['E', 'm'])

def hamiltonian(Q, 
                eigvals, eigvecs, contour, 
                num_particles=2, verbose=False):
    sep_M = n_n.gen_matrix(eigvecs, contour, Q)
    
    # tuples of k-indexes: (3, 6), (4, 4)
    mb_E = list(combinations_with_replacement(Q.E, num_particles))
    order = len(mb_E)
    
    def H_func(i, j):
            return coupled_H_elem(mb_E[i], mb_E[j], eigvals, sep_M, Q)
    return matrix_from_function(H_func, order)

def coupled_H_elem(E_bra, E_ket, eigvals, sep_M, Q):
    E1, E2 = E_ket
    E1_prim, E2_prim = E_bra
    mod = 1
    if E1 == E2:
        mod = mod/2
    if E1_prim == E2_prim:
        mod = mod/2
    
    one_body = 0
    if E_bra == E_ket:
        one_body += eigvals[E1] + eigvals[E2]
    
    two_body = sep_M[E1, E1_prim] * sep_M[E2, E2_prim] + sep_M[E1, E2_prim] * sep_M[E2, E1_prim]
    two_body = n_n.V0 * two_body
    
    return mod * (one_body + two_body)

def coupled_H_elem_prev(E_bra, E_ket, eigvals, Q):
    E1, E2 = E_bra
    E1_prim, E2_prim = E_ket
    
    result = 0
    for m in Q.m:
        bra = FermionState([SP(m, E1), SP(Q.M - m, E2)])
        if bra.sign is 0: continue
        
        for m_prim in Q.m:
            ket = FermionState([SP(m_prim, E1_prim), 
                                SP(Q.M - m_prim, E2_prim)])
            if ket.sign is 0: continue
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
            one_body += eigvals[sp_state.E]
    
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







"""
TESTS
"""
    
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

