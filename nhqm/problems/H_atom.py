from __future__ import division
import scipy as sp

name = "H_atom"
mass = 1
HO_omega = 1
eV_factor = 27.2113
energy_units = "Hartrees"

@sp.vectorize
def V(r, l, s):
    return - 1. / r

potential = V