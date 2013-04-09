from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations
import mom_space as mom
from collections import defaultdict
import many_body as mb


def gen_states(s, num_particles):
    return mb.gen_states(s,num_particles)

def get_sp_smart(num_states, num_particles, res_dict):
    mb_states= gen_states(num_states, num_particles)
    sp_states= gen_states(num_states, num_particles=1)
        
    for bra_index, bra in enumerate(mb_states):
        for alfa in bra:
            alpha=set([alfa])
            for beta in (sp_states): 
                """we'd like to generate beta from the set of possible sp_states minus the set of states already present in the (bra -alpha), can this be done through clever set manipulations?"""
                if ((bra - alpha) - beta) == (bra - alpha):
                    #can't add existing fermion 
                    ket = (bra - alpha) | beta
                    ket_index = mb_states.index(ket)
                    res = (alpha, beta)
                    res_dict[bra_index, ket_index][0].append(res)
    return  
        
"""def translate_ket_to_index( state_vector, ket):
    return state_vector.index(ket)        """
        
def get_2p_smart(num_states, num_particles, res_dict):
    mb_states = gen_states(num_states, num_particles)
    twop_states = gen_states(num_states, num_particles = 2)
    
    
    for bra_index, bra in enumerate(mb_states):
        for alphabeta_interim in combinations(bra, 2):
            alphabeta = set(alphabeta_interim)
            for gammadelta in twop_states: 
                """ we'd like t generate gammadelta from all twop states minus any two-particle-state formed from two of the particles already present in the (bra -alphabeta), can this be done through clever set manipulations?"""
                ket = (bra - alphabeta)
                if ket - gammadelta == ket:
                    #can't add existing fermions   
                    ket = ket | gammadelta
                    ket_index = mb_states.index(ket)
                    res = (alphabeta, gammadelta) 
                    res_dict[(bra_index,ket_index)][1].append( res )
    return

def get_hamilton_dict(ns,np):
    def default_factory():
        return ([],[])
    h_dict = defaultdict( default_factory )
    get_2p_smart(ns,np,h_dict)
    get_sp_smart(ns,np,h_dict)
    return h_dict
    
                 
    
if __name__ == '__main__':
        

    ns =5
    np =3
    
    res = get_hamilton_dict(ns,np) 
    
    sps = gen_states(ns,np)
    a,b,c,d,e,f = 0,0,2,3,2,4
    print sps[a], sps[b], "::", res[(a,b)][0]
    print sps[c], sps[d], "::", res[(c,d)][0]
    print sps[e], sps[f], "::", res[(e,f)][0]
    print ""
    print "2p:"
    print ""
    print sps[a], sps[b], "::", res[(a,b)][1]
    print sps[c], sps[d], "::", res[(c,d)][1]
    print sps[e], sps[f], "::", res[(e,f)][1]
    