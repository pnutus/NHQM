from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom

V0 = -100 #MeV
beta = 1 # fm^-2

def interaction(a, b, c, d):
    return sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c]

def gen_matrix(eigvecs, contour, Q):
    order = len(eigvecs)
    global sep_M
    sep_M = sp.empty((order, order), complex)
    for i, bra in enumerate(eigvecs):
        for j, ket in enumerate(eigvecs):
            
            ### is this were we want to calculate one element per gpu kernel?
            sep_M[i, j] = separable_elem(bra, ket, contour, Q)
    sep_M *= sp.sqrt(V0)
    return sep_M

def separable_elem(bra, ket, contour, Q):
    zip_contour = zip(*contour)
    result = 0
    #no weight in the contour array?
    for n, (k, w) in enumerate(zip_contour):
        inner_sum = 0
        #Complex conjugate the bra? Or NHQM trickery?
        
        ### should this innerproduct sum be evaluated on only one kernel or would it significantly improve performance if this ~50 element dot/inner product was in itself parallelized?
        for n_prim, (k_prim, w_prim) in enumerate(zip_contour):
            inner_sum += w_prim * ket[n_prim] * V_sep(k, k_prim, Q.l, Q.j)   
        result += w * bra[n] * inner_sum 
    return result

def V_sep(k, k_prim, l, j):
    args = (k, k_prim, potential, l, j)
    
    ### should this integral be the run on one kernel (albeit in paralell with other integrals like this one) or split up between many kernels? seeing as it' only about 30 points would there be much performance to gain or is it of lower importance?
    integral, _ = fixed_quad(mom.integrand, 0, 10, n = 20, args=args)
    return 2 / sp.pi * k_prim**2 * integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
# def mom.integrand(r, k, k_prim, V, l, j):
#     return r**2 * j_l(l, k*r) * j_l(l, k_prim*r) * V(r, l, j)