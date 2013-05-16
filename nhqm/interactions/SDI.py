from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function
from nhqm.wigner import wigner3j

name = "Delta"

V0 = -1670 # MeV
r0 = 2 # fm

def gen_interaction(eigvecs, Q, sp_basis, contour=None):
    if sp_basis == osc:
        R = [osc.gen_wavefunction(eigvec)(r0) 
                for eigvec in eigvecs.T]
    elif sp_basis == mom:
        R = [mom.gen_wavefunction(eigvec, contour, Q)(r0)
                for eigvec in eigvecs.T]
    hatjs = (2*Q.j + 1)**2
    wigners = wigner3j(Q.j, Q.j, Q.J, 0.5, 0.5, -1)**2 \
            + wigner3j(Q.j, Q.j, Q.J, 0.5, -0.5, 0)**2
                  
    def interaction(a, b, c, d):
        Kabcd = V0*r0**2*R[a]*R[b]*R[c]*R[d] / (16*sp.pi)
        return Kabcd * hatjs * wigners
    return interaction