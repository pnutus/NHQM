from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, combinations_with_replacement
from nhqm.bases.fermion_state import FermionState, gen_mb_states
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function

name = "Uncoupled"

def hamiltonian(Q, eigvals, eigvecs, sep_M, num_particles=2):
    # Generate many-body states
    mb_states = gen_mb_states(Q, num_particles)
    # Only take states with m1 + m2 = M
    def total_M(mb_state):
        return Q.M == sum(state.m for state in mb_state)
    mb_states = filter(total_M, mb_states)
    
    # Generate hamiltonian
    matrix_size = len(mb_states)
    def H_func(i, j):
        return H_elem(mb_states[i], mb_states[j], eigvals, sep_M)
    return matrix_from_function(H_func, matrix_size, symmetric=True)
            
def H_elem(bra, ket, eigvals, sep_M):
    # H_1, one-body interaction
    one_body = 0
    if bra == ket:
        one_body = sum(eigvals[sp.E] for sp in bra)
    
    # H_2, two-body interaction
    a, b = bra
    c, d = ket
    two_body = 0
    if a.m == c.m and b.m == d.m:
        two_body += sep_M[a.E, c.E]*sep_M[b.E, d.E]
    if a.m == d.m and b.m == c.m:
        two_body -= sep_M[a.E, d.E]*sep_M[b.E, c.E]
        
    return one_body + two_body
