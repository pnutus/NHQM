from imports import *

def V(r, l, j):
    """Simplest possible potential."""
    return 0

class CentralProblem:
    """Describes a quantum mechanical problem with spherical symmetry."""
    
    def __init__(self, name):
        self.name = name
        self.potential = V # Spherically symmetric potential
        self.mass = 1 # Mass in units of the problem
        self.eV_factor = 27.3 # Scaling factor to get to eV
    

        