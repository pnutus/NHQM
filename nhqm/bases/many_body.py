from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, permutations, imap

# THIS IS FOR FERMIONS

def gen_states(num_sp_states, num_particles=1):
    if num_sp_states < num_particles:
        raise ValueError(
        "There cannot be more particles than states"
        )
    return imap(set, combinations(range(num_sp_states), num_particles))

def H_elem(bra, ket, sp_H, eigvecs):
    one_body = sum(sp_H[a, b] 
                    for (a, b) in one_body_indexes(bra, ket))
    two_body = sum(n_n_interaction(a, b, c, d, eigvecs)
                    for (a, b, c, d) in two_body_indexes(bra, ket))
    return one_body + two_body
    
def n_n_interaction(a, b, c, d, eigvecs):
    """
    Some kinda integration I guess
    """

def one_body_indexes(bra, ket):
    result = []
    for b in ket:
        for a in bra:
            new_ket = ket - set([b]) | set([a])
            if new_ket == bra:
                result.append( (a, b) )
    return result

def two_body_indexes(bra, ket):
    result = []
    for annihilated in combinations(ket, 2):
        for created in combinations(bra, 2):
            new_ket = (ket
                        .difference(annihilated)
                        .union(created))
            if new_ket == bra:
                result.append( created + annihilated )
    return result

# Write (regression) tests! Also verify!

# states = set([0, 3, 1]), set([0, 3, 4])
# print one_body_indexes(*states)
# print two_body_indexes(*states)
# print gen_states(30*4, 2)