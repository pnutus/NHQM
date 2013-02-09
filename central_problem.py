from imports import *
import scipy.constants as const

def V(r, l, j):
    """Simplest possible potential."""
    return 0

class CentralProblem:
    """Describes a quantum mechanical problem with spherical symmetry."""
    
    def __init__(self, name):
        self.name = name
        self.potential = V # Spherically symmetric potential
        self.mass = 1 # Mass in units of the problem
        # Scaling factor to get to eV
        self.eV_factor = const.physical_constants["atomic unit of energy"][0] \
                         / const.eV # hartree