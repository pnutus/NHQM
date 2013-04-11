from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations
import mom_space as mom
from collections import defaultdict
import many_body as mb
from fermion_state import FermionState


def gen_states(s, num_particles):
    return mb.gen_states(s,num_particles)
    
def gen_separable_matrix(eigvecs, contour):
    return mb.gen_separable_matrix(eigvecs, contour) 
    
def n_n_interaction(sep_M, a, b, c, d):
    return sep_M[a, c]*sep_M[a, d] - sep_M[a, b]*sep_M[d, c]         


def hamiltonian(eigvals, eigvecs, contour, num_particles=2):
    num_sp_states = len(eigvecs)
    mb_states = gen_states(num_sp_states, num_particles)
    order = len(mb_states)
    sep_M = gen_separable_matrix(eigvecs, contour)
    H = sp.zeros( (order,order), complex )
    hamilton_dict = get_hamilton_dict(num_sp_states, num_particles)
    
    for key, values in hamilton_dict.iteritems():
        twop_interactions = values
        "two body; n-n interaction:"  
        for alphabeta, gammadelta in twop_interactions:
            ab = list(alphabeta)
            gd = list(gammadelta)
            
            #print mb_states[key[0]], mb_states[key[1]], ab, gd
            temp = n_n_interaction(sep_M, int(ab[0]), int(ab[1]),\
             int(gd[0]), int(gd[1]) )
            H[key] += temp
    return H
        
def find_n_n_interactions(num_states, num_particles, res_dict):
    mb_states = gen_states(num_states, num_particles)
    twop_states = gen_states(num_states, num_particles = 2)
    
    
    for bra_index, bra in enumerate(mb_states):
        for alphabeta in combinations(bra, 2):
            
            
            for gammadelta in twop_states: 
                """ we'd like t generate gammadelta from all 
                twop states minus any two-particle-state formed 
                from two of the particles already present in 
                the (bra -alphabeta), can this be done through 
                clever set manipulations?"""
                ket = bra.annihilate(alphabeta).create(gammadelta)
                if ket.sign != 0:
                    #can't add existing fermions 
                    ket_index = mb_states.index(FermionState(ket))
                    res = (FermionState(alphabeta), gammadelta, ket.sign) 
                    res_dict[(bra_index,ket_index)].append( res )
    return

def get_hamilton_dict(num_states, num_particles):

    """h_dict contains the value [] between two bra 
    and ket indices that do not permitt a fock space 
    contribution, idealy the key shouldn't even be 
    represented in the dict if neither sp nor 2p 
    permitts a contribution. and if only one does so 
    then no [] value from the other shold be saved"""    
        
    h_dict = defaultdict( list )
    find_n_n_interactions(num_states, num_particles, h_dict)
    return h_dict

            
    
if __name__ == '__main__':
    print "poor grammar makes me [sic]"