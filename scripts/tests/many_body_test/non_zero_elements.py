from imports import *
from nhqm.bases import many_body as mb
from nhqm.bases import naive_many_body as nmb
from nhqm.bases import many_body_spill as mbs
from nhqm.bases.fermion_state import FermionState

import scipy as sp
from scipy import linalg
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import itertools
            
import unittest

def mb_get_interactions(num_state, num_part):
    result = []
    mb_states = mb.gen_states(num_state, num_part)
    for bra in mb_states:
        for ket in mb_states:
            for res in mb.two_body_indexes(bra, ket):
                result.append((mb_states.index(bra), mb_states.index(ket), res))
    return result            
   
class RedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print "starting set up"
        cls.num_s = 10
        cls.num_p = 4

        #cls.naive = nmb.get_two_particle_combinations(
        #                cls.num_s, cls.num_p)               
                        
        cls.medium = mb_get_interactions(
                        cls.num_s, cls.num_p)
                        
        cls.smart = mbs.get_hamilton_dict(cls.num_s, cls.num_p)
        
        cls.mb_s = mb.gen_states(cls.num_s, cls.num_p)
        
        #print "medium:", cls.medium
        #print " smart:", cls.smart               
        
    def testLength(self):
        
        len_smart = 0
        for key,values in RedTests.smart.iteritems():
            for element in values:
                len_smart += 1
    
        
        print "medium length:", len(RedTests.medium)
        print " smart length:", len_smart
        print ""
        
        
        
        self.assertEquals(len(RedTests.medium), len_smart)
         
    def testElements(self):
        check = 0
        
        for i in xrange(len(RedTests.medium)):
            
            medium_key = (RedTests.medium[i][0], RedTests.medium[i][1])
            #print cls.medium[i]
            medium_ab = RedTests.medium[i][2][0:2]
            medium_cd = RedTests.medium[i][2][2:4]
            #print medium_key, medium_abcd_list
            
            for key,values in RedTests.smart.iteritems():
                for element in values:
                    if key == medium_key: 
                        smart_abcd_list = element[0].states + element[1].states
                        smart_abcd_list.append(element[2])
                        if smart_abcd_list == list(RedTests.medium[i][2]):
                            check += 1
                        
                        
        self.assertEquals(check, len(RedTests.medium))                
        
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

