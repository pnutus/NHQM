
from __future__ import division
from imports import *
from nhqm.bases import many_body as mb
from nhqm.bases import many_body_spill as mbs
import scipy as sp
from scipy import linalg
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import itertools
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
from time import time
            
import unittest


def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))


problem = He5.problem   
order = 3*1
l = 1
j = 1.5
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 2.5

contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
#points, _ = contour
zip_contour = zip(*contour)
args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)
idx = res_index(eigvecs)
res_E = eigvals[idx]

Ham = mbs.hamiltonian(eigvals, eigvecs, contour, num_particles=2)
print Ham