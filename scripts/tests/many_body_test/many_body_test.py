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
                print res
                result.append((bra, ket, res))
    return result            
   
class RedTests(unittest.TestCase):
    
    def setUp(self):
        print "starting set up"
        self.num_s = 3
        self.num_p = 2

        self.naive = nmb.get_two_particle_combinations(
                        self.num_s, self.num_p)               
                        
        self.medium = mb_get_interactions(
                        self.num_s, self.num_p)
                        
        self.smart = mbs.get_hamilton_dict(self.num_s, self.num_p)
        
        self.mb_s = mb.gen_states(self.num_s, self.num_p)
        
    """def testMedium2PLength(self):
        N = len(self.naive) - len(self.medium)
        self.assertEquals(N, 0 )  
        
        
    def testMedium2PElements(self):
        check =0
        length = len(self.medium)
        for i in xrange(length):
            if self.medium[i] == self.naive[i]:
                check +=1 
                   
        self.assertEquals(length, check)"""                   
        
        
        
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

