from __future__ import division
from imports import *
from itertools import combinations_with_replacement
from collections import namedtuple
from nhqm.problems import He5
from nhqm.QM_helpers import energies, symmetric, hermitian
from nhqm.bases.gen_contour import triangle_contour, gauss_contour
from nhqm.bases import (mom_space            as mom, 
                        harm_osc             as osc, 
                        mb_coupled           as coupled, 
                        mb_uncoupled         as uncoupled,
                        two_body_interaction as n_n)
from nhqm.plot_helpers import *
import numpy as np

#finds bound state and resonance
problem = He5 
basis_size = 15*3
n_n.V0 = -1140
n_n.r0 = 1
peak_x = 0.5
peak_y = 0.5
k_max = 12
problem.V0 = -47.
real_contour = False



#problem = He5 
#basis_size = 25*3
#n_n.V0 = -5810
#n_n.r0 = 0.5
#peak_x = 0.6
#peak_y = 0.5
#k_max = 20
#problem.V0 = -47.
#real_contour = False

#problem = He5 
#basis_size = 10*3
#n_n.V0 = -98
#n_n.r0 = 2
#peak_x = .5
#peak_y = .5
#k_max = 5
#problem.V0 = -47.
#real_contour = False


contour = triangle_contour(peak_x, peak_y, k_max, basis_size/3)
if real_contour:
    contour = gauss_contour([0, k_max], basis_size)
points, _ = contour
QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))

def main():
    print ("Real contour" if real_contour else "Complex Contour")
    print "n-n V0:", n_n.V0
    print "Basis size:", basis_size
    print "k_max:", k_max
    solve_3b(mom, coupled)
    # solve_3b(mom, uncoupled)
    # solve_3b(osc, coupled)
    # solve_3b(osc, uncoupled)

def solve_3b(sp_basis, mb_scheme):
    global eigvals
    eigvals, eigvecs, sep_M = solve_2b(sp_basis)
    mb_H = mb_scheme.hamiltonian(Q, eigvals, eigvecs, 
                                 sep_M, num_particles = 2)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    
    print mb_scheme.name, sp_basis.name, 
    print "lowest energy:" , mb_eigvals[0]
    plot_shit(eigvals, mb_eigvals, mb_eigvecs)
    plt.show()

def solve_2b(basis):
    if basis == osc:
        eigvals, eigvecs = osc.solution(basis_size, problem, Q)
    else:
        eigvals, eigvecs = mom.solution(contour, problem, Q)
    sep_M = n_n.gen_matrix(eigvecs, Q, basis, contour)
    return eigvals, eigvecs, sep_M
    
def plot_shit(eigvals, mb_eigvals, mb_eigvecs):
    def res_index():
        result = 0
        global_min = 1
        for i in range(len(mb_eigvecs)):
            temp_max = max(abs(mb_eigvecs[:,i]))
            if temp_max < global_min:
                global_min = temp_max
                result = i
        return result
    res = res_index()
    
    plt.figure(1)
    plt.clf()
    #plot_poles(mb_eigvals, problem.mass)
    ks = sp.sqrt(2*problem.mass*mb_eigvals)
    plt.plot(sp.real(ks), sp.imag(ks), 'ko')
    
    k = sp.zeros((len(points), len(points)),complex)
    for i, E1 in enumerate(eigvals):
        for j, E2 in enumerate(eigvals):
            k[i,j]=sp.sqrt(2*problem.mass*(E1 + E2))
    plt.plot(sp.real(k),sp.imag(k),'or', markersize=4)
    plt.axis([0, 2.1 * 1.42 * peak_x, - 1.1 * 1.42 * peak_y, 0])
    
    plt.figure(2)
    plt.clf()
    #mb_E = list(combinations_with_replacement(Q.E, 2))
    #Es = sp.zeros(len(mb_E), complex)
    #Es_2 = sp.zeros((len(mb_E),2), complex)
    #for i, (E_1, E_2) in enumerate(mb_E):
    #    Es[i] = eigvals[E_1]+eigvals[E_2]
    #    Es_2[i,:] = (eigvals[E_1], eigvals[E_2])
    #    #Es[i] = max(abs(eigvals[E_1]), abs(eigvals[E_2]))
    #args = np.argsort(mb_eigvecs[:,res])
    #plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,0])),'b')
    #plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,1])),'g')
    #plt.figure(3)
    #plt.clf()
    #plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,0])),abs(mb_eigvecs[:,res]),'b*')
    #plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,1])),abs(mb_eigvecs[:,res]),'g*')
    plt.figure(4)
    plt.clf()
    plt.plot(abs(mb_eigvecs[:,res]))
    print 'energy =', mb_eigvals[res], 'MeV'
    print 'momentum =', sp.sqrt(2*problem.mass*(mb_eigvals[res]))

if __name__ == '__main__':
    main()
