
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
            
import unittest

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

class RedTests(unittest.TestCase):
    
    def setUp(self):
        print "starting set up"

        problem = He5.problem   
        order = 5
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
        
        self.mb_H, self.sep  = mb.hamiltonian(H, eigvecs, zip_contour, num_particles=2)
        self.mbs_H, self.seps = mbs.hamiltonian(H, eigvecs, zip_contour, num_particles=2)
        
        self.matrix_row = self.mb_H.shape[0]
        self.matrix_col = self.mb_H.shape[1]
        
        num_sp_states = len(eigvecs)
        self.mb_s = mb.gen_states(num_sp_states,2)
        
        print self.sep == self.seps
        
    def testMatrixEquiv(self):
        print "matrix elements"
        N=0
        
        for i in xrange(self.matrix_row):
            for j in xrange(self.matrix_col):
                if abs(self.mb_H[i,j] - self.mbs_H[i,j]) < 0.1:
                    N+=1
                else:
                    
                    print self.mb_s[i], self.mb_s[j], ":::"
                    print "mb: ", self.mb_H[i,j]
                    print "mbs:", self.mbs_H[i,j]  
                    print "" 
        self.assertEquals(N, self.matrix_row*self.matrix_col )
        

if __name__ == '__main__':
    unittest.main()        