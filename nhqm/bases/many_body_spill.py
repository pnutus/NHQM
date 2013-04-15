from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations, product
from collections import namedtuple
import mom_space as mom
from collections import defaultdict
import many_body as mb
from fermion_state import FermionState


def state_index(states, new_state):
    for i, state in enumerate(states):
        if new_state == state:
            raise ValueError
        if new_state < state: #lt?
            return i
    return len(states)

def gen_mb_states(quantum_numbers, num_particles=2):
    """
    Generates many-body states given a list of 
    quantum numbers and their possible values.
    Like so [('m', [-1.5, -0.5, 0.5, 1.5]),
             ('k', [0, 1, 2, 3, 4, 5, 6])]
    """
    names, values = zip(*quantum_numbers)
    sp = namedtuple('sp', names)
    sp_states = [sp(*tup) for tup in product(*values)]
    return map(FermionState, combinations(sp_states, num_particles))


    
def gen_separable_matrix(eigvecs, contour):
    return mb.gen_separable_matrix(eigvecs, contour) 
    
def n_n_interaction(sep_M, a, b, c, d):
    return sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c]

def hamiltonian(eigvals, eigvecs, contour, num_particles=2):
    

    q_nums = [('m', [-1.5, -0.5, 0.5, 1.5]),
              ('k', range( len(contour) ))]
              
    mb_states = gen_mb_states(q_nums, 2)
    
    num_sp_states = len(eigvecs)
    
    order = len(mb_states)
    sep_M = gen_separable_matrix(eigvecs, contour)
    H = sp.zeros( (order,order), complex )
    hamilton_dict = get_hamilton_dict(num_sp_states, num_particles)
    
    
    
    #one body energy:    
    for i, value in enumerate(mb_states):
        a, b = value
        
        H[i,i] = eigvals[a] + eigvals[b]
    
    
    #two body; n-n interaction:
    for key, twop_interactions in hamilton_dict.iteritems():
        for alphabeta, gammadelta, sign in twop_interactions:
            ab = list(alphabeta)
            gd = list(gammadelta)
            
            #print mb_states[key[0]], mb_states[key[1]], ab, gd
            res = n_n_interaction(sep_M, int(ab[0]), int(ab[1]),\
             int(gd[0]), int(gd[1]) )
            H[key] += res
    return H
        
def find_n_n_interactions(num_states, num_particles, res_dict):
    
    q_nums = [('m', [-1.5, -0.5, 0.5, 1.5]),
              ('k', range(num_particles))]
    mb_states = gen_mb_states(q_nums, 2)   #already calculated in hamiltonian...
    
    q2_nums = [('m', [-1.5, -0.5, 0.5, 1.5]),
              ('k', range(2))]
    twop_states = gen_mb_states(q2_nums, 2)
    print mb_states
    
    for bra_index, bra in enumerate(mb_states):
        for alphabeta in combinations(bra, 2):
            
            
            for gammadelta in twop_states: 
                print bra
                print alphabeta
                print gammadelta
                
                ket = bra.annihilate(alphabeta).create(gammadelta)
                
                print ket
                
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