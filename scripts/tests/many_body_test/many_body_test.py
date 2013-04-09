from imports import *
from nhqm.bases import many_body as mb
from nhqm.bases import naive_many_mody as nmb
from scipy.sparse import lil_matrix
import scipy as sp
from scipy import linalg
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import itertools
            
import unittest
   
class RedTests(unittest.TestCase):
    
    def setUp(self):
        print "starting set up"
        self.num_s = 7
        self.num_p = 3
        self.naive_sp = nmb.get_single_particle_combinations(
                        self.num_s, self.num_p)
        #self.smart_sp = mb.get_single_particle_combinations(
                        #self.num_s, self.num_p)
        self.smart_sp = mb.get_sp_smart(self.num_s, self.num_p)

        self.naive_2p = nmb.get_two_particle_combinations(
                        self.num_s, self.num_p)
        #self.smart_2p = mb.get_two_particle_combinations(
        #                self.num_s, self.num_p)
        self.smart_2p = mb.get_2p_smart(self.num_s, self.num_p)
        
    def testSPLength(self):
        print "testing single particle length"
        N = len(self.naive_sp) - len(self.smart_sp)
        print "n",len(self.naive_sp)
        print "s", len(self.smart_sp)
        self.assertEquals(N, 0 )

    def test2PLength(self):
        print "testing two particle length"
        N = len(self.naive_2p) - len(self.smart_2p)
        print "s2", len(self.smart_2p)
        print "n2", len(self.naive_2p)
        self.assertEquals(N, 0 )
    
    def testSPElements(self):
        print "correlating single particle sorting order"
        length = min( len(self.naive_sp), len( self.smart_sp ) )
        check_2 = 0
        
        self.smart_sp.sort()
        self.naive_sp.sort()
        
        for i in xrange(length):
            for j in xrange(length):
                if self.smart_sp[i] == self.naive_sp[j]:
                    check_2 = check_2 + 1
        
        #for i in xrange(length):
        #    if not ( self.smart_sp[i] == self.naive_sp[i] ):
        #        num_errors = num_errors +1
        self.assertEquals(length, check_2)

    def test2PElements(self):
        print "correlating two particle sorting order"
        length = min(len(self.naive_2p), len( self.smart_2p ) )
        res, ref = sp.zeros(length), sp.zeros(length)
        check = 0
        
        for i in xrange(length):
            for j in xrange(length):
                if self.naive_2p[i] == self.smart_2p[j]:
                    check = check +1
        self.assertEquals(check, length)

        
class ComplexTests(unittest.TestCase):
    pass
    
class WeightTests(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()
    
    
    
    
    
    
    
    
    
    
    
    
    
"""problem = He5.problem   
order = 1*3
l = 0
j = 0.5
problem.V0 = -47.
peak_x = 0.17
peak_y = 0.07
k_max = 2.5
num_states = 5
num_particles = 2

contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
points, _ = contour

args = (problem, l, j)
H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
eigvals, eigvecs = calc.energies(H)

eigvecs, _ = linalg.eig(H)
states = mb.gen_states(order, num_particles)
mb_H = mb.hamiltonian(states, H, eigvecs)


###

naivecomb = nmb.get_naive_combinations_single_particle(num_states, num_particles)
comb = mb.get_single_combinations(num_states, num_particles)

num_d = 10
p_d = 3        

naivedoublecomb = nmb.get_naive_combinations_two_particle_interaction(num_d, p_d)  

doublecomb=mb.get_two_particle_combinations(num_d, p_d)

b = len(naivedoublecomb) 
a = len(doublecomb)
print b - a

for i in xrange( min( len(naivedoublecomb), len(doublecomb) ) ):
    if naivedoublecomb[i] != doublecomb[i]:
        print naivedoublecomb[i], doublecomb[i]
        
"""

