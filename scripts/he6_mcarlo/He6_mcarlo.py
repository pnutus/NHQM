from __future__ import division
# from imports import *
from itertools import combinations_with_replacement
from collections import namedtuple
from nhqm.problems import He5
from nhqm.QM_helpers import energies, symmetric, hermitian
from nhqm.bases.contours import triangle_contour_explicit, gauss_contour
from nhqm.bases import mom_space as mom, harm_osc as osc
from nhqm.interactions import gaussian, SDI
from nhqm.mb_schemes import coupled, uncoupled
from nhqm.plot_helpers import *
import scipy
import numpy
import random
#import resonances_2body as res2b

sp_basis = mom
mb_scheme = coupled
two_body = SDI

# SDI trickery
problem = He5 
problem.V0 = -47.
basis_size = 45
points_on_triangle = basis_size*2/3
gaussian.V0 = -1140
gaussian.r0 = 1
SDI.V0 = 2800
SDI.r0 = 2
peak_x = 0.2
peak_y = 0.1
k_max = 2.5
complex_contour = True 
J = 0
js = [1.5]

if complex_contour:
    contour = triangle_contour_explicit(peak_x, peak_y, k_max, 
                    points_on_triangle, basis_size - points_on_triangle)
else:
    contour = gauss_contour([0, k_max], basis_size)
points, _ = contour

print mb_scheme.name, sp_basis.name
print "\tBasis size:", basis_size
if sp_basis == mom:
    print "\t",("Complex contour" if complex_contour else "Real contour")
    print "\t\tk_max:", k_max
print two_body.name, "interaction"
print "\tV0:", two_body.V0
print "\tr0:", two_body.r0


# Construct sp states.
spQ = namedtuple('qnums', 'l j')
SP = namedtuple('sp', "id l j E eigvec contour basis")
sp_states = []
for j in js:
    # We can use different contours here.
    Q = spQ(l=1, j=j)
    eigvals, eigvecs = sp_basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        sp_states.append(SP(id=len(sp_states), l=1, j=j, E=eigvals[i], 
                eigvec=eigvecs[:,i], contour=contour, basis=sp_basis))
# for i,sp in enumerate(sp_states):
    # e = scipy.sqrt(2*problem.mass*sp.E)
    # print i, ":", e
    # plt.plot(scipy.real(e),scipy.imag(e),'*')
# plt.show()

# Helper function.
def unshared_copy(inList):
    if isinstance(inList, list):
        return list( map(unshared_copy, inList) )
    return inList

# Helper function.
def get_rand(seg,res,states):
    rand = -1
    while rand == -1 or rand == res or rand in states:
        rand = random.sample(range(int(seg*basis_size/3),int((seg+1)*basis_size/3)),1)[0]
    return rand

n_rounds=100
resonance_i=14
n_p_segs=1
results=[]
for i in range(n_rounds):
    states = []
    for i in range(3*n_p_segs):
        states.append(sp_states[get_rand(i%3,resonance_i,states)])
    states.append(sp_states[resonance_i])
    for i,sp in enumerate(states):
        states[i] = SP(i,l=sp.l,j=sp.j,E=sp.E,eigvec=sp.eigvec,contour=sp.contour,basis=sp.basis)
    interaction = two_body.gen_interaction(states, J, problem=problem)
    mb_H = coupled.hamiltonian(states, interaction, J)
    mb_eigvals, mb_eigvecs = energies(mb_H)
    results.append(mb_eigvals)

ks=[]
for e in results:
    ks.append(scipy.sqrt(2*problem.mass*e))
plt.plot(scipy.real(ks),scipy.imag(ks),'*')
plt.show()