from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function


V0 = -1000 # MeV
r0 = 2 # fm

def gen_interaction(eigvecs, Q, sp_basis, contour=None):
    
    if sp_basis == osc:
        R = [osc.gen_wavefunction(eigvec)(r0) 
                for eigvec in eigvecs.T]
    elif sp_basis == mom:
        R = [mom.gen_wavefunction(eigvec, contour, Q)(r0)
                for eigvec in eigvecs.T]
    def interaction(a, b, c, d):
        return V0 * r0**2 * R[a.E] * R[b.E] * R[c.E] * R[d.E]
    return interaction