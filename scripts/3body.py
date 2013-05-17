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
#import resonances_2body as res2b

sp_basis = mom
mb_scheme = coupled
two_body = SDI

# SDI trickery
problem = He5 
problem.V0 = -47.
basis_size = 15
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


# Construct sp states
spQ = namedtuple('qnums', 'l j')
SP = namedtuple('sp', "id l j E eigvec contour basis")
sp_states = []
for j in js:
    # We can use different contours here
    Q = spQ(l=1, j=j)
    eigvals, eigvecs = sp_basis.solution(contour, problem, Q)
    for i in xrange(len(eigvals)):
        sp_states.append(SP(id=len(sp_states), l=1, j=j, E=eigvals[i], 
                eigvec=eigvecs[:,i], contour=contour, basis=sp_basis))

# for state in sp_states:
#     print state.id, state.j, state.E

interaction = two_body.gen_interaction(sp_states, J, problem=problem)
mb_H = coupled.hamiltonian(sp_states, interaction, J)
mb_eigvals, mb_eigvecs = energies(mb_H)

print "Lowest energy:" , mb_eigvals[0]

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
    args = sp.argsort(mb_eigvecs[:,res])
    plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,0])),'b')
    plt.plot(abs(sp.sqrt(2*problem.mass*Es_2[args,1])),'g')
    plt.figure(3)
    plt.clf()
    plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,0])),abs(mb_eigvecs[:,res]),'b*')
    plt.plot(sp.sqrt(abs(2*problem.mass*Es_2[:,1])),abs(mb_eigvecs[:,res]),'g*')
    plt.figure(4)
    plt.clf()
    plt.plot(abs(mb_eigvecs[:,res]))
    print 'Maybe resonance energy =', mb_eigvals[res], 'MeV'
    print 'Maybe resonance momentum =', sp.sqrt(2*problem.mass*(mb_eigvals[res]))

# plot_shit(eigvals, mb_eigvals, mb_eigvecs)
# plt.show()



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
