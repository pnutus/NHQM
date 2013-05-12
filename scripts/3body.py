from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import triangle_contour
from nhqm.bases import many_body as mb
from nhqm.bases import two_body_interaction as n_n
from nhqm.plot_helpers import *
from itertools import combinations_with_replacement

problem = He5 
order = 10*3

problem.V0 = -47.
peak_x = 0.17
peak_y = 0.2
k_max = 3

n_n.V0 = -5000


contour = triangle_contour(peak_x, peak_y, k_max, order/3)
points, _= contour

Q = mom.QNums(l=1, j=1.5, k=range(len(contour[0])))
      
H = mom.hamiltonian(contour, problem, Q)
eigvals, eigvecs = energies(H)

Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
             m=[-1.5, -0.5, 0.5, 1.5], 
             E=range(len(eigvals)))
             
mb_H = mb.hamiltonian(Q, eigvals, eigvecs, contour, num_particles = 2, verbose=True)
mb_eigvals, mb_eigvecs = energies(mb_H)
print "lowest energy:", mb_eigvals[0]

plt.figure(1)
plt.clf()
plot_poles(mb_eigvals, problem.mass)
#plt.plot(sp.real(mb_eigvals), sp.imag(mb_eigvals), 'ko')

k = sp.zeros((len(points), len(points)),complex)
for i, ks1 in enumerate(points):
    for j, ks2 in enumerate(points):
        k[i,j]=sp.sqrt(ks1 ** 2 + ks2 ** 2)
plt.plot(sp.real(k),sp.imag(k),'og')
        
#plt.show()

"""
ONE-BODY TESTS
"""
    
import unittest
class RedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        
        problem = He5 
        order = 1*3
        problem.V0 = -47.
        peak_x = 0.17
        peak_y = 0.2
        k_max = 4

        n_n.V0 = 0

        contour = triangle_contour(peak_x, peak_y, k_max, order/3)
        points, _= contour

        Q = mom.QNums(l=1, j=1.5, k=range(len(contour[0])))
      
        H = mom.hamiltonian(contour, problem, Q)
        eigvals, eigvecs = energies(H)

        Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
                     m=[-1.5, -0.5, 0.5, 1.5], 
                     E=range(len(eigvals)))
        num_particles = 2     
        mb_H = mb.hamiltonian(Q, eigvals, eigvecs, contour, num_particles, verbose=True)
        mb_eigvals, mb_eigvecs = energies(mb_H)
        
        cls.hamilton_test =mb_H
        cls.eigvals = eigvals 
        cls.order = order
        cls.mb_E = list(combinations_with_replacement(Q.E, num_particles))
        cls.mod = True
    
    def setUp(self):
        print ""
        
        
    def test00(self):
        a = 0
        b = 0
        c = self.mb_E.index((a,b))
        if self.mod and a == b: 
            self.assertEquals(4 * self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )
        else:
            self.assertEquals(self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )    
        
    def test01(self):
        a = 0
        b = 1
        c = self.mb_E.index((a,b))        
        if self.mod and a == b: 
            self.assertEquals(4 * self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )
        else:
            self.assertEquals(self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )    
        
    def test11(self):
        a = 1
        b = 1
        c = self.mb_E.index((a,b))
        if self.mod and a == b: 
            self.assertEquals(4 * self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )
        else:
            self.assertEquals(self.hamilton_test[c,c], self.eigvals[a] + self.eigvals[b] )    
        
#if __name__ == '__main__':
    #print "kaptenkvant 4 lyfe"
    unittest.main()    
