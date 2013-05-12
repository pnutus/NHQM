from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom, \
                       harm_osc as osc
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import triangle_contour
from nhqm.bases import  many_body as mom_coupled, \
                        many_body_uncoupled as mom_uncoupled, \
                        many_body_uncoupled_HO as osc_uncoupled
from nhqm.bases import two_body_interaction as n_n
from nhqm.plot_helpers import *
from itertools import combinations_with_replacement

problem = He5 
problem.V0 = -47.

n_n.V0 = -500
order = 25*3

Q = mom_coupled.QNums(l=1, j=1.5, J=0, M=0, 
                      m=[-1.5, -0.5, 0.5, 1.5], 
                      E=range(order))

peak_x = 0.17
peak_y = 0.2
k_max = 4
contour = triangle_contour(peak_x, peak_y, k_max, order/3)
points, _= contour

def main():
    # run_mom_coupled()
    # run_mom_uncoupled()
    run_osc_uncoupled()

def run_mom_coupled():
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)

    mb_H = mom_coupled.hamiltonian(Q, eigvals, eigvecs, 
                                contour, num_particles = 2)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    print "Coupled MomSpace lowest energy:", mb_eigvals[0]

def run_mom_uncoupled():
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    
    mb_H = mom_uncoupled.hamiltonian(Q, eigvals, eigvecs, 
                                contour, num_particles = 2)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    print "Uncoupled MomSpace lowest energy:", mb_eigvals[0]

def run_osc_uncoupled():
    H = osc.hamiltonian(order, problem, Q)
    eigvals, eigvecs = energies(H)
    
    mb_H = osc_uncoupled.hamiltonian(Q, eigvals, eigvecs, 
                                    num_particles = 2)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    print "Uncoupled HarmOsc lowest energy:", mb_eigvals[0]
    
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
        
    plt.show()

if __name__ == '__main__':
    main()