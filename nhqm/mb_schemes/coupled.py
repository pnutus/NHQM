from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations_with_replacement
from collections import namedtuple, Iterable
from nhqm.QM_helpers import matrix_from_function, energies
 
name = "Coupled"

def solution(sp_states, interaction, J):
    H = hamiltonian(sp_states, interaction)
    return energies(H)

def hamiltonian(sp_states, interaction, J):
    # tuples of sp states (sp(E=1, j=1.5), sp(E=7, j=1.5))
    mb_states = list(combinations_with_replacement(sp_states, 2))
    matrix_size = len(mb_states)
    
    def H_func(i, j):
        return H_elem(mb_states[i], mb_states[j], interaction, J)
    return matrix_from_function(H_func, matrix_size)

def H_elem(bra, ket, interaction, J):
    a, b = ket
    c, d = bra
    
    one_body = 0
    if bra == ket:
        one_body = sum(sp.E for sp in bra) * (1 + (a == b)*(-1)**J)
        
    two_body = interaction(a, b, c, d)
    
    return (one_body + two_body) * N(a, b, J) * N(c, d, J)

def N(a, b, J):
    return sp.sqrt(1 + (a == b)*(-1)**J)/(1 + (a == b))