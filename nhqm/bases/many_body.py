from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from itertools import combinations
import mom_space as mom
from fermion_state import FermionState

sqrtV0 = sp.sqrt(-2250) #sqrt(MeV)
beta = 1 # fm^-2
# THIS IS FOR FERMIONS

def hamiltonian(eigvals, eigvecs, contour, num_particles=2, verbose=False):
    if verbose: print "Generating many-body states"
    num_sp_states = len(eigvecs)
    mb_states = gen_states(num_sp_states, num_particles)
    order = len(mb_states)
    if verbose: print "Generating matrix of separable interactions"
    sep_M = gen_separable_matrix(eigvecs, contour)
    if verbose: print "Generating many-body H"
    H = sp.empty((order, order), complex)
    for i, bra in enumerate(mb_states):
        for j, ket in enumerate(mb_states):
            H[i,j] = H_elem(bra, ket, eigvals, sep_M)
    return H, sep_M
    
def gen_states(num_sp_states, num_particles=2):
    if num_sp_states < num_particles:
        raise ValueError("There cannot be more particles than states")
    return map(FermionState, combinations(range(num_sp_states), num_particles))  

def H_elem(bra, ket, eigvals, sep_M):
    if bra == ket:
        a, b = bra
        one_body = eigvals[a] + eigvals[b]
    else:
        one_body = 0
        
    def n_n_interaction(a, b, c, d):
        return sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c]
    two_body = sum(sign * n_n_interaction(a, b, c, d)
                    for (a, b, c, d, sign) in two_body_indexes(bra, ket)) 
                                  
    return one_body + two_body

def gen_separable_matrix(eigvecs, contour):
    order = len(eigvecs)
    M = sp.empty( (order, order), complex )
    for i, bra in enumerate(eigvecs):
        for j, ket in enumerate(eigvecs):
            M[i, j] = separable_elem(bra, ket, contour)
    return M

def separable_elem(bra, ket, contour, l=0, j=0.5):
    zip_contour = zip(*contour)
    result = 0
    #no weight in the contour array?
    for n, (k, w) in enumerate(zip_contour):
        inner_sum = 0
        #Complex conjugate the bra? Or NHQM trickery?
        for n_prim, (k_prim, w_prim) in enumerate(zip_contour):
            inner_sum += w_prim * ket[n_prim] * V_sep(k, k_prim, l, j)
        result += w * bra[n] * inner_sum 
    return result

def V_sep(k, k_prim, l, j):
    args = (k, k_prim, potential, l, j)
    integral, _ = fixed_quad(mom.integrand, 0, 10, n = 20, args=args)
    return 2 / sp.pi * k_prim**2 * integral
    
def potential(r, l, j):
    return sqrtV0*sp.exp(- beta * r**2)

def two_body_indexes(bra, ket):
    result = []
    if len(set(bra) - set(ket)) > 2:
        return []
    for annihilated in combinations(ket, 2):
        for created in combinations(bra, 2):
            new_ket = (ket
                        .annihilate(annihilated)
                        .create(created))
            if new_ket.states == bra.states:
                sign = new_ket.sign
                result.append(created + annihilated + (sign,))
    return result
    
if __name__ == '__main__':
    print "kaptenkvant 4 lyfe"
    print two_body_indexes(FermionState([1,2,3]), FermionState([1,3,4]))
