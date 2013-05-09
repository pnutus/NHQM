from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, combinations_with_replacement
from collections import namedtuple
from nhqm.bases.fermion_state import FermionState, gen_mb_states
import nhqm.bases.two_body_interaction as n_n
from nhqm.QM_helpers import clebsch_gordan, matrix_from_function
            


# THIS IS FOR FERMIONS

QNums = namedtuple('qnums', 'l J M j E m')

def hamiltonian(Q, eigvals, eigvecs, contour, 
                num_particles=2, verbose=False):
    sep_M = n_n.gen_matrix(eigvecs, contour, Q)
    
    def total_M(mb_state):
        return Q.M == sum(state.m for state in mb_state)
    # tuples of k-indexes: (3, 6), (4, 4)
    mb_states = gen_mb_states(Q, num_particles)
    mb_states = filter(total_M, mb_states)
    order = len(mb_states)
    
    def H_func(i, j):
            return H_elem(mb_states[i], mb_states[j], eigvals, sep_M)
    return matrix_from_function(H_func, order)
            
def H_elem(bra, ket, eigvals, sep_M):
    # H_1, one-body interaction
    one_b =0
    twob_1 =0
    twob_2 =0
    if bra == ket:
        one_b = sum(eigvals[sp.E] for sp in bra)
  
    
    # H_2, two-body interaction
    bra1, bra2 = bra
    ket1, ket2 = ket
    if bra1.m == ket1.m and bra2.m == ket2.m:
        twob_1 =  n_n.V0 * sep_M[bra1.E, ket1.E] * sep_M[bra2.E, ket2.E]
    if bra1.m == ket2.m and bra2.m == ket1.m:
        twob_2 = - n_n.V0 * sep_M[bra1.E, ket2.E] * sep_M[bra2.E, ket1.E]
    
    if not True and one_b == twob_2 == twob_1 == 0:
        print bra, ket
        print one_b, twob_2, twob_1    
        print "::::"
    return one_b + twob_2 + twob_1
