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

def get_naive_combinations_single_particle(num_states, num_particles=1):

    states=gen_states(num_states, num_particles)
    states2=gen_states(num_states, num_particles)
    single_state=gen_states(num_states, num_particles=1)
        
    combs = []        

    for i, a in enumerate(states):
        for j, b in enumerate(states2):
                #testa alfa och beta och se om det blir nagot
                for beta in b:
                    for alpha in single_state:
                        if (a == (b - set([beta]) | alpha )):
                            combs.append( \
                            (a, b, alpha, set([beta])) )
                        
    return combs                       
            
def get_naive_combinations_two_particle_interaction(num_states, num_particles=2):

    states=gen_states(num_states, num_particles)
    states2=gen_states(num_states, num_particles)
    double_state=gen_states(num_states, num_particles=2)
        
    combs = []        

    for i, a in enumerate(states):
        for j, b in enumerate(states2):
                #testa alfa1,2 och beta1,2 och se om det blir nagot
                gamma_delta_permut = combinations(b,2)
                for gamma_delta in gamma_delta_permut:
                    for alpha_beta in double_state:
                        if (a == (b - set(gamma_delta) | alpha_beta )):
                            combs.append( \
                            (a, b, alpha_beta, set(gamma_delta) ) )
                        
    return combs                        
            
            
import unittest
   
class RedTests(unittest.TestCase):
    
    def setUp(self):
        self.num_s = 10
        self.num_p = 3
        self.naive_sp = mb.get_single_particle_combinations(
                        self.num_s, self.num_p)
        self.smart_sp = nmb.get_single_particle_combinations(
                        self.num_s, self.num_p)

        self.naive_2p = mb.get_two_particle_combinations(
                        self.num_s, self.num_p)
        self.smart_2p = nmb.get_two_particle_combinations(
                        self.num_s, self.num_p)
    
    def testSPLength(self):
        N = len(self.naive_sp) - len(self.smart_sp)
        self.assertEquals(N, 0 )

    def testSPLength(self):
        N = len(self.naive_2p) - len(self.smart_2p)
        self.assertEquals(N, 0 )
    
    def testSPElements(self):
        length = min( len(self.naive_sp), len( self.smart_sp ) )
        num_errors = 0
        
        for i in xrange(length):
            if not ( self.smart_sp[i] == self.naive_sp[i] ):
                num_errors = num_errors +1
        self.assertEquals(num_errors, 0)

    def test2PElements(self):
        length = min(len(self.naive_2p), len( self.smart_2p ) )
        res, ref = sp.zeros(length), sp.zeros(length)
        num_errors = 0
        
        for i in xrange(length):
            if not( self.smart_2p[i] == self.naive_2p[i] ):
                num_errors = num_errors +1
        self.assertEquals(num_errors, 0)

        
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


