from __future__ import division
from imports import *
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

problem = He5 
basis_size = 10*3
n_n.V0 = -1800
peak_x = 0.17
peak_y = 0.2
k_max = 4
problem.V0 = -47.

contour = triangle_contour(peak_x, peak_y, k_max, basis_size/3)
#contour = gauss_contour([0, k_max], basis_size)
points, _ = contour
QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))

def main():
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
    plot_shit(mb_eigvals)
    plt.show()

def solve_2b(basis):
    if basis == osc:
        eigvals, eigvecs = osc.solution(basis_size, problem, Q)
    else:
        eigvals, eigvecs = mom.solution(contour, problem, Q)
    sep_M = n_n.gen_matrix(eigvecs, Q, basis, contour)
    return eigvals, eigvecs, sep_M
    
def plot_shit(mb_eigvals):
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
if __name__ == '__main__':
    main()
