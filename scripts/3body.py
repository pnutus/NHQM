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
problem.V0 = -47.

n_n.V0 = -1200
basis_size = 7*3

QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(basis_size))

peak_x = 0.17
peak_y = 0.2
k_max = 5
contour = triangle_contour(peak_x, peak_y, k_max, basis_size/3)
contour = gauss_contour([0, k_max], basis_size)
points, _ = contour

def main():
    solve_3b(mom, coupled)
    solve_3b(mom, uncoupled)
    solve_3b(osc, coupled)
    solve_3b(osc, uncoupled)

def solve_3b(sp_basis, mb_scheme):
    eigvals, eigvecs, sep_M = solve_2b(sp_basis)
    
    mb_H = mb_scheme.hamiltonian(Q, eigvals, eigvecs, 
                                 sep_M, num_particles = 2)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    
    print mb_scheme.name, sp_basis.name, 
    print "lowest energy:" , mb_eigvals[0]

def solve_2b(basis):
    if basis == osc:
        H = osc.hamiltonian(basis_size, problem, Q)
    else:
        H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    sep_M = n_n.gen_matrix(eigvecs, Q, basis, contour)
    return eigvals, eigvecs, sep_M
    
def plot_shit(mb_eigvals):
    plt.figure(1)
    plt.clf()
    plot_poles(mb_eigvals, problem.mass)
    #plt.plot(sp.real(mb_eigvals), sp.imag(mb_eigvals), 'ko')

    k = sp.zeros((len(points), len(points)),complex)
    for i, ks1 in enumerate(points):
        for j, ks2 in enumerate(points):
            k[i,j]=sp.sqrt(ks1 ** 2 + ks2 ** 2)
    plt.plot(sp.real(k),sp.imag(k),'og')

if __name__ == '__main__':
    main()