from __future__ import division
from imports import *
from itertools import combinations_with_replacement
from collections import namedtuple
from nhqm.problems import He5
from nhqm.QM_helpers import energies, symmetric, hermitian
from nhqm.bases.contours import triangle_contour_explicit, gauss_contour
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.interactions import gaussian, SDI
from nhqm.mb_schemes import coupled, uncoupled
from nhqm.plot_helpers import *
import numpy as np
#import resonances_2body as res2b

# SDI trickery
problem = He5 
problem.V0 = -47.
basis_size = 75
points_on_triangle = 10
gaussian.V0 = -1140
gaussian.r0 = 1
SDI.V0 = -1670
SDI.r0 = 2
peak_x = 0.3
peak_y = 0.3
k_max = 30
complex_contour = False

if complex_contour:
    contour = triangle_contour_explicit(peak_x, peak_y, k_max, 
                                        points_on_triangle, basis_size - points_on_triangle)
else:
    contour = gauss_contour([0, k_max], basis_size)
points, _ = contour
QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))

def main():
    solve_3b(mom, coupled, SDI)

def solve_3b(sp_basis, mb_scheme, two_body):
    print mb_scheme.name, sp_basis.name
    print "\tBasis size:", basis_size
    if sp_basis == mom:
        print "\t",("Complex contour" if complex_contour else "Real contour")
        print "\t\tk_max:", k_max
    print two_body.name, "interaction"
    print "\tV0:", two_body.V0
    print "\tr0:", two_body.r0
    eigvals, eigvecs = solve_2b(sp_basis)
    interaction = two_body.gen_interaction(eigvecs, Q, sp_basis, contour)
    mb_H = mb_scheme.hamiltonian(Q, eigvals, eigvecs, 
                                 interaction, num_particles = 2)
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
    # p12_en, p12_vec = res2b.export_resonance_p12(contour)
    # sp.append(eigvals,p12_en)
    # sp.append(eigvecs,p12_vec)
    return eigvals, eigvecs
    
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
    #plt.axis([0, 2.1 * 1.42 * peak_x, - 1.1 * 1.42 * peak_y, 0])
    
    plt.figure(2)
    plt.clf()
    mb_E = list(combinations_with_replacement(Q.E, 2))
    Es = sp.zeros(len(mb_E), complex)
    Es_2 = sp.zeros((len(mb_E),2), complex)
    for i, (E_1, E_2) in enumerate(mb_E):
        Es[i] = eigvals[E_1]+eigvals[E_2]
        Es_2[i,:] = (eigvals[E_1], eigvals[E_2])
        #Es[i] = max(abs(eigvals[E_1]), abs(eigvals[E_2]))
    args = np.argsort(mb_eigvecs[:,res])
    plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,0])),'b')
    plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,1])),'g')
    plt.figure(3)
    plt.clf()
    plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,0])),abs(mb_eigvecs[:,res]),'b*')
    plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,1])),abs(mb_eigvecs[:,res]),'g*')
    plt.figure(4)
    plt.clf()
    plt.plot(abs(mb_eigvecs[:,res]))
    print 'energy =', mb_eigvals[res], 'MeV'
    print 'momentum =', sp.sqrt(2*problem.mass*(mb_eigvals[res]))

if __name__ == '__main__':
    main()




#finds bound state and resonance
# problem = He5 
# basis_size = 15*3
# gaussian.V0 = -1140
# gaussian.r0 = 1
# peak_x = 0.5
# peak_y = 0.5
# k_max = 12
# k_max = 20
# problem.V0 = -47.
# complex_contour = True

# problem = He5 
# basis_size = 25*3
# gaussian.V0 = -5810
# gaussian.r0 = 0.5
# peak_x = 0.6
# peak_y = 0.5
# k_max = 20
# problem.V0 = -47.
# complex_contour = True
  
# problem = He5 
# basis_size = 10*3
# gaussian.V0 = -98
# gaussian.r0 = 2
# peak_x = .5
# peak_y = .5
# k_max = 5
# problem.V0 = -47.
# complex_contour = True
