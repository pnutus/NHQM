from imports import *
from nhqm.bases import many_body as mb
from nhqm.bases import many_body_spill as mbs
import scipy as sp
from scipy import linalg
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import itertools
from __future__ import division
from plotting import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
            
import unittest

class RedTests(unittest.TestCase):
    
    def setUp(self):
        print "starting set up"

        problem = He5.problem   
        order = 50
        l = 1
        j = 1.5
        problem.V0 = -47.
        peak_x = 0.17
        peak_y = 0.07
        k_max = 2.5

        contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
        points, _ = contour
        args = (problem, l, j)
        H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
        eigvals, eigvecs = calc.energies(H)
        idx = res_index(eigvecs)
        res_E = eigvals[idx]
        
        mb_H = mb.hamiltonian(H, eigvecs, contour, num_particles=2)
        mbs_H = mbs.hamiltonian(H, eigvecs, contour, num_particles=2)
        
    def testSPLength(self):
        print "matrix elements"
        N=0
        
        for i in xrange(dim(mb_h)[0]):
            for j in xrange(dim(mb_h)[1]):
                if abs(mb_H[i,j] - mbs_H[i,j]) > 0.000001:
                    N+=1
        self.assertEquals(N, 0 )