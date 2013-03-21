from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations
import mom_space as mom

# THIS IS FOR FERMIONS

def hamiltonian(sp_H, eigvecs, contour, num_particles=2):
    num_sp_states = len(eigvecs)
    mb_states = gen_states(num_sp_states, num_particles)
    order = len(mb_states)
    sep_M = gen_separable_matrix(eigvecs, contour)
    H = sp.empty( (order,order), complex )
    for i, bra in enumerate( mb_states ):
        for j, ket in enumerate( mb_states ):
            H[i,j] = H_elem(bra, ket, sp_H, sep_M)
    return H  
    
def gen_states(num_sp_states, num_particles=2):
    if num_sp_states < num_particles:
        raise ValueError("There cannot be more particles than states")
    return map(set, combinations(range(num_sp_states), num_particles))

def H_elem(bra, ket, sp_H, sep_M):
    V0 = 1
    def n_n_interaction(a, b, c, d):
        return sep_M[a, c]*sep_M[a, d] - sep_M[a, b]*sep_M[d, c]
    one_body = sum(sp_H[a, b] 
                    for (a, b) in one_body_indexes(bra, ket))
    two_body = sum(n_n_interaction(a, b, c, d)
                    for (a, b, c, d) in two_body_indexes(bra, ket))
    return one_body + two_body

def gen_separable_matrix(eigvecs, contour):
    order = len(eigvecs)
    M = sp.empty( (order, order), complex )
    for i, bra in enumerate(eigvecs):
        for j, ket in enumerate(eigvecs):
            M[i, j] = separable_elem(bra, ket, contour)
    return M

def separable_elem(bra, ket, contour, l=0, j=0.5):
    result = 0
    for n, (k, w) in enumerate(contour):
        inner_sum = 0
        for n_prim, (k_prim, w_prim) in enumerate(contour):
            inner_sum += w_prim * ket[n_prim] * V_sep(k, k_prim, l, j)
        result += w * bra[n] * inner_sum 
    return result

def V_sep(k, k_prim, l, j):
    args = (k, k_prim, potential, l, j)
    integral, _ = fixed.quad(mom.integrand, 0, 10, n = 20, args=args)
    return 2 / sp.pi * k_prim**2 * integral
    
def potential(r, l, j):
    beta = 1
    sqrtV0 = 1
    return sqrtV0*sp.exp(- beta * r**2)
    
def one_body_indexes(bra, ket, verbose=False):
    result = []
    for b in ket:
        for a in bra:
            new_ket = ket - set([b]) | set([a])
            if new_ket == bra:
                if verbose:
                    result.append( (bra, ket, set([a]), set([b]) ) )
                else:    
                    result.append( (a, b) )
    return result

def two_body_indexes(bra, ket,verbose = False):
    result = []
    for annihilated in combinations(ket, 2):
        for created in combinations(bra, 2):
            new_ket = (ket
                        .difference(annihilated)
                        .union(created))
            if new_ket == bra:
                if verbose:
                    result.append( 
                            (bra, ket, set(created), set(annihilated) ) )
                else:
                    result.append( created + annihilated )
    return result
    

def get_single_particle_combinations(num_states, num_particles = 1):
    mb_states = gen_states(num_states, num_particles)
    result = []
    for i, bra in enumerate( mb_states ):
        for j, ket in enumerate( mb_states ):
            temp =  one_body_indexes(bra,ket,verbose=True) 
            
            if temp != []:
                for k, tup in enumerate(temp):
                 result.append( tup )
    return result        
            
def get_two_particle_combinations(num_states, num_particles=2):
    mb_states = gen_states(num_states, num_particles)
    result = []
    for i, bra in enumerate(mb_states):
        for j, ket in enumerate(mb_states):
            temp = two_body_indexes(bra,ket,verbose=True)
            
            if temp != []:
                for k, tup in enumerate(temp):
                 result.append( tup )
    
    return result    
    

