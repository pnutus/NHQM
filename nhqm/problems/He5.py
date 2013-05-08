from __future__ import division
import scipy as sp

name = "He5"
mass = 0.019272
eV_factor = 1e6
units = "MeV"
HO_omega = 20

V0 = -47 # MeV
Vso = -7.5 # MeV
r0 = 2. # fm
d = .65 # fm

def V(r, l, j):
    f = 1/(1 + sp.e**((r - r0)/d))
    spin_orbit = .5*(j*(j + 1) - l*(l + 1) - .75)
    return f * (V0 - 4*Vso*spin_orbit*(f - 1) / (d * r))

potential = V