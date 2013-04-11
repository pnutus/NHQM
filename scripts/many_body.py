from __future__ import division
from imports import *
from plotting import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
from nhqm.bases import many_body as mb
from nhqm.bases.fermion_state import FermionState

mb.sqrtV0 = sp.sqrt(-2139) #sqrt(MeV)

problem = He5.problem   
order = 5*3
l = 1
j = 1.5
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 4

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
points, _ = contour
args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)

H, sep_M = mb.hamiltonian(eigvals, eigvecs, contour, num_particles = 2, verbose=True)
eigvals, eigvecs = calc.energies(H)
print "lowest energy:", eigvals[0]

# plt.pcolor(H)
# plt.gca().invert_yaxis()
# plt.show()