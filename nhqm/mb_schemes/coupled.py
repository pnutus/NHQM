from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations_with_replacement
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function
 
name = "Coupled"

def hamiltonian(Q, eigvals, eigvecs, interaction, num_particles=2):
    
    # tuples of k-indexes: (3, 6), (4, 4)
    mb_E = list(combinations_with_replacement(Q.E, num_particles))
    order = len(mb_E)
    
    def H_func(i, j):
        return coupled_H_elem(mb_E[i], mb_E[j], eigvals, interaction, Q)
    return matrix_from_function(H_func, order)

def coupled_H_elem(E_bra, E_ket, eigvals, interaction, Q):
    # H_1, one-body interaction
    one_body = 0
    a, b = E_ket
    c, d = E_bra
    if E_bra == E_ket:
        one_body = sum(eigvals[E] for E in E_bra) * (1 + (a == b)*(-1)**Q.J)
    
    # H_2, two-body interaction
    two_body = (interaction(a, b, c, d) 
              + interaction(a, b, d, c) * (-1)**Q.J )
    
    return (one_body  + two_body) * N(a, b, Q) * N(c, d, Q)

def N(a, b, Q):
    return sp.sqrt(1 + (a == b)*(-1)**Q.J)/(1 + (a == b))
