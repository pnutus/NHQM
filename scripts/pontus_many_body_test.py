from __future__ import division
from imports import *
from plotting import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
from nhqm.bases import many_body as mb

problem = He5.problem   
order = 10*3
l = 1
j = 1.5
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 4

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))
        
# Old for comparison:
#contour = calc.naive_triangle_contour(0, 0, k_max, order)
contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
points, _ = contour
args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)

H = mb.hamiltonian(eigvals, eigvecs, contour)
eigvals, eigvecs = calc.energies(H)
print eigvals