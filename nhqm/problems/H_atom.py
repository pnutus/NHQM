from imports import *
from central_problem import CentralProblem

@sp.vectorize
def V(r, l, s):
    return - 1. / r
    
problem = CentralProblem("Hydrogen atom")
problem.potential = V
problem.mass = 1
problem.HO_omega = 1
problem.eV_factor = 27.2113
problem.units = "Hartrees"
