from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, permutations

# THIS IS FOR FERMIONS

def gen_states(num_sp_states, num_particles=1):
    if num_sp_states < num_particles:
        raise ValueError(
        "There cannot be more particles than states"
        )
    return map(set, combinations(range(num_sp_states), num_particles))
    
def gen_state_mask(num_sp_states, num_particles =1):
    
    states = gen_states(num_sp_states, num_particles)
    res =  sp.zeros( (len(states),num_sp_states) )
    for i, state in enumerate(states):
        for j in state:
            res[i,j]=1
    return res               

def hamiltonian(mb_states, sp_H, eigvecs):
    order = len(mb_states)
    H = sp.empty( (order,order), complex )
    for i, bra in enumerate( mb_states ):
        for j, ket in enumerate( mb_states ):
            H[i,j] = H_elem(bra, ket, sp_H, eigvecs)
            
    return H            

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
    return 0;
    
def one_body_indexes(bra, ket):
    result = []
    for b in ket:
        for a in bra:
            new_ket = ket - set([b]) | set([a])
            if new_ket == bra:
                result.append( (a, b) )
    return result

def one_body_comb(bra, ket):
    result = []
    for b in ket:
        for a in bra:
            new_ket = ket - set([b]) | set([a])
            if new_ket == bra:
                result.append( (bra, ket) )
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

def get_combinations(num_states, num_particles = 1):
    mb_states = gen_states(num_states, num_particles)
    result = []
    for i, bra in enumerate( mb_states ):
        for j, ket in enumerate( mb_states ):
            temp =  one_body_comb(bra,ket) 
            
            if temp != []:
                for k, tup in enumerate(temp):
                 result.append( tup )
    return result        
            
# Write (regression) tests! Also verify!

if __name__ == '__main__':
    mask = gen_state_mask(5,2)
    a = mask[0,:]
    b = mask[1,:]
#     states = set([0,1,3]), set([0,1,5])
#     print one_body_indexes(*states)
#     print two_body_indexes(*states)
#     print gen_states(30, 2)
