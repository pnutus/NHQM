from imports import *
import central_problem import CentralProblem

@sp.vectorize
def V(r, l, s):
    return - 1. / r
    
problem = CentralProblem()
problem.potential = V
problem.mass = 1
problem.eV_factor = 13.6