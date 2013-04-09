from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations
import mom_space as mom
from collections import defaultdict
import many_body as mb


def gen_states(s, num_particles):
    return mb.gen_states(s,num_particles)
    
def gen_separable_matrix(eigvecs, contour):
    return mb.gen_separable_matrix(eigvecs, contour)    


def hamiltonian(sp_H, eigvecs, contour, num_particles=2):
    num_sp_states = len(eigvecs)
    mb_states = gen_states(num_sp_states, num_particles)
    order = len(mb_states)
    sep_M = gen_separable_matrix(eigvecs, contour)
    H = sp.empty( (order,order), complex )
    hamilton_dict = get_hamilton_dict(num_sp_states, num_particles)
    
    for key in hamilton_dict.keys():
        "single particle energy:"
        for alfa_beta_tup in hamilton_dict[key][0]:
            alpha = alfa_beta_tup[0]
            beta = alfa_beta_tup[1]
            
            H[key] += sp_H[int(alfa),int(beta)]
        "two body; n-n interaction:"  
        for alphabeta_gammadelta_tup in hamilton_dict[key][1]:
            ab = alphabeta_gammadelta_tup[0]
            gd= alphabeta_gammadelta_tup[1]
            
            H[key] += mb.n_n_interaction( int(ab[0]), int(ab[1]),\
             int(gd[0]), int(gd[1]) )
    
    return H


def get_sp_smart(num_states, num_particles, res_dict):
    mb_states= gen_states(num_states, num_particles)
    sp_states= gen_states(num_states, num_particles=1)
    
    #res_n =0    
    
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
                    
                    #print res_n
                    #res_n += 1
                    #print res
                    
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
    
    """h_dict contains the value [] between two bra and ket indices that do not permitt a fock space contribution, idealy the key shouldn't even be represented in the dict if neither sp nor 2p permitts a contribution. and if only one does so then no [] value from the other shold be saved"""    
        
    h_dict = defaultdict( default_factory )
    get_2p_smart(ns,np,h_dict)
    get_sp_smart(ns,np,h_dict)
    return h_dict

def get_sp_dict(ns,np):
    def default_factory():
        return ([],[])
    h_dict = defaultdict( default_factory )
    get_sp_smart(ns,np,h_dict)
    return h_dict  
      
def get_2p_dict(ns,np):
    def default_factory():
        return ([],[])
    h_dict = defaultdict( default_factory )
    get_2p_smart(ns,np,h_dict)
    return h_dict                   
    
if __name__ == '__main__':
        

    ns =5
    np =3
    
    res = get_hamilton_dict(ns,np) 
    
    mb_s = gen_states(ns,np)
    a,b,c,d,e,f = 0,0,2,3,2,4
    print mb_s[a], mb_s[b], "::", res[(a,b)][0]
    print mb_s[c], mb_s[d], "::", res[(c,d)][0]
    print mb_s[e], mb_s[f], "::", res[(e,f)][0]
    print ""
    print "2p:"
    print ""
    print mb_s[a], mb_s[b], "::", res[(a,b)][1]
    print mb_s[c], mb_s[d], "::", res[(c,d)][1]
    print mb_s[e], mb_s[f], "::", res[(e,f)][1]
    
    print ""
    print""
    print ""
    print ":::::::::::::::::::::::::::"
    
    #for key in res:
    #    print mb_s[key[0]], mb_s[key[1]], "::", res[key]