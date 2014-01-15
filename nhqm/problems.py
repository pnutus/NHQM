from __future__ import division
import scipy as sp

class Helium5:
    # hbar = 1, energy in MeV, length in fm
    name = "Helium5"
    mass = 0.019272
    eV_factor = 1e6 # = 1 MeV
    energy_units = "MeV"
    HO_omega = 20 # close to optimal, improves harm osc basis expansion
    
    def __init__( self
                , V0 = -47 # MeV
                , Vso = -7.5 # MeV
                , r0 = 2. # fm
                , d = .65 # fm
                ):
        self.V0  = V0
        self.Vso = Vso
        self.r0  = r0
        self.d   = d
       
    def potential(self, r, l, j):
        f = 1/(1 + sp.e**((r - self.r0)/self.d))
        spin_orbit = .5*(j*(j + 1) - l*(l + 1) - .75)
        return f * (self.V0 - 4*self.Vso*spin_orbit*(f - 1) / (self.d * r))

class HydrogenAtom:
    # atomic units
    name = "Hydrogen Atom"
    mass = 1
    HO_omega = 1 # close to optimal, improves harm osc basis expansion
    eV_factor = 27.2113 # = 1 Hartree
    energy_units = "Hartrees"
    
    def __init__(self):
        pass
    
    @staticmethod
    def potential(r, l, j):
        return - 1. / r