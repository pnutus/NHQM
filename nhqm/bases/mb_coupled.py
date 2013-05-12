from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations_with_replacement
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function
            
name = "Coupled"

def hamiltonian(Q, eigvals, eigvecs, sep_M, num_particles=2):
    
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
        mod /= sp.sqrt(2)
    if E1_prim == E2_prim:
        mod /= sp.sqrt(2)
    
    one_body = 0
    if E_bra == E_ket:
        one_body = eigvals[E1] + eigvals[E2]
    
    two_body =  (sep_M[E1, E1_prim] * sep_M[E2, E2_prim] 
           + 1 * sep_M[E1, E2_prim] * sep_M[E2, E1_prim])
    
    return mod * (one_body + two_body)
   
        