from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
import mom_space as mom

V0 = -1e9 #MeV
beta = 1 # fm^-2

def interaction(a, b, c, d):
    return V0*(sep_M[a, c]*sep_M[b, d] - sep_M[a, d]*sep_M[b, c])

def gen_matrix(eigvecs, contour, Q, verbose=False):
    order = len(eigvecs)
    global sep_M
    sep_M = sp.empty((order, order), complex)
    for i, bra in enumerate(eigvecs):
        for j, ket in enumerate(eigvecs):
            if verbose: print "{}%".format((i*order + j)/order*100)
            sep_M[i, j] = separable_elem(bra, ket, contour, Q)
    sep_M *= 2 / sp.pi
    return sep_M

def separable_elem(bra, ket, contour, Q):
    zip_contour = zip(*contour)
    result = 0
    #no weight in the contour array?
    for n, (k, w) in enumerate(zip_contour):
        inner_sum = 0
        #Complex conjugate the bra? Or NHQM trickery?
        for n_prim, (k_prim, w_prim) in enumerate(zip_contour):
            inner_sum += w_prim * ket[n_prim] * k_prim**2 \
                            * V_sep(k, k_prim, Q.l, Q.j)
        result += w * bra[n] * inner_sum * k**2
    return result

def V_sep(k, k_prim, l, j):
    args = (k, k_prim, potential, l, j)
    integral, _ = fixed_quad(mom.integrand, 0, 10, n = 20, args=args)
    return integral

@sp.vectorize
def potential(r, l, j):
    return sp.exp(- beta * r**2)
    
