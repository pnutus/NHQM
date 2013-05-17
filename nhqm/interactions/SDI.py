from __future__ import division
import scipy
from scipy.integrate import fixed_quad
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.QM_helpers import matrix_from_function
from nhqm.wigner import wigner3j
from itertools import product

#Only for uncoupled!

name = "Delta"

V0 = 1670 # MeV
r0 = 2 # fm

def gen_interaction(sp_states, J, problem=None):
    R = [sp.basis.gen_wavefunction(sp.eigvec, sp, 
                        problem=problem, contour=sp.contour)(r0)
                        for sp in sp_states]
    
    wigner = {}
    js = set(sp.j for sp in sp_states)
    for j1, j2 in product(js, js):
        wigner[(j1, j2)] = wigner3j(j1, j2, J, 0.5, -0.5, 0)
    
    def hatj(sp):
        return 2*sp.j + 1
    
    # from suhonen
    def interaction(a, b, c, d):
        if J % 2 == 1:
            return 0
        Kabcd = - V0 * r0**2 * R[a.id]*R[b.id]*R[c.id]*R[d.id] / (16*scipy.pi)
        wigners = wigner[(a.j, b.j)] * wigner[(c.j, d.j)]
        hatjs = scipy.sqrt(hatj(a)*hatj(b)*hatj(c)*hatj(d))
        return - Kabcd * (-1)**(b.j + d.j) * 4 * hatjs * wigners
    return interaction
