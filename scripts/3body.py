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
problem.V0 = -47.05
problem.Vso = -7.04
basis_size = 40
points_on_triangle = 30
gaussian.V0 = -7140
gaussian.r0 = 1
SDI.V0 = 998
SDI.r0 = 2
peak_x = 0.3
peak_y = 0.3
k_max = 2.5
complex_contour = True
J = 2
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

he6_spstates=[]
for sp_state in sp_states:
    if sp.real(sp.sqrt(2*problem.mass*sp_state.E)) < 0.7:
        he6_spstates.append(sp_state)
he6_spstates = sp_states

interaction = two_body.gen_interaction(sp_states, J, problem=problem)
mb_H = coupled.hamiltonian(he6_spstates, interaction, J)
mb_eigvals, mb_eigvecs = energies(mb_H)

print "Lowest energy:" , mb_eigvals[0]

def plot_shit(eigvals, mb_eigvals, mb_eigvecs):
    
    def find_index(vector, value):
        for x,y in enumerate(vector):
            if y==value:
                return x
    
    def res_index(eigvecs):
        maxes = map(max, abs(eigvecs.T))
        return maxes.index(min(maxes))
    
    sp_res = res_index(eigvecs)
    mb_states = list(combinations_with_replacement(he6_spstates, 2))
    Es = sp.zeros((len(mb_states),2), complex)
    for i, (sp_state_1, sp_state_2) in enumerate(mb_states):
        Es[i,:] = (sp_state_1.E, sp_state_2.E)
        if sp_state_1.E == eigvals[sp_res] and sp_state_2.E == eigvals[sp_res]:
            dbl_res = i
    
    temp_sorted = sp.sort(abs(mb_eigvecs[dbl_res,:]))
    res_1 = find_index(abs(mb_eigvecs[dbl_res,:]), temp_sorted[-1])
    res_2 = find_index(abs(mb_eigvecs[dbl_res,:]), temp_sorted[-2])
    print 'Maybe resonance energy =', mb_eigvals[res_1], 'MeV'
    print 'Maybe resonance momentum =', sp.sqrt(2*problem.mass*(mb_eigvals[res_1]))
    plt.figure(1)
    plt.clf()
    #plot_poles(mb_eigvals, problem.mass)
    ks = sp.sqrt(2*problem.mass*mb_eigvals)
    plt.plot(sp.real(ks), sp.imag(ks), 'ko')
    plt.plot(sp.real(ks[res_1]), sp.imag(ks[res_1]), 'go', markersize=10)
    plt.plot(sp.real(ks[res_2]), sp.imag(ks[res_2]), 'bo', markersize=10)
    
    print len(mb_states)
    import os
    script_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
    dataoladata = sp.array([sp.real(ks), sp.imag(ks)])
    sp.savetxt(script_dir + "he6_momenta"+ str(res_1)+".data", dataoladata.T)

    k = sp.zeros((len(points), len(points)),complex)
    for i, E1 in enumerate(eigvals):
        for j, E2 in enumerate(eigvals):
            k[i,j]=sp.sqrt(2*problem.mass*(E1 + E2))
    plt.plot(sp.real(k),sp.imag(k),'or', markersize=4)
    #plt.axis([0, 0.3, -0.02, 0])
    
    plt.figure(2)
    plt.clf()

    
        #Es[i] = max(abs(eigvals[E_1]), abs(eigvals[E_2]))
    args = sp.argsort(mb_eigvecs[:,res_1])
    plt.plot(abs(sp.sqrt(2*problem.mass*Es[args,0])),'b')
    plt.plot(abs(sp.sqrt(2*problem.mass*Es[args,1])),'g')
    plt.figure(3)
    plt.clf()
    plt.plot(sp.sqrt(abs(2*problem.mass*Es[:,0])),abs(mb_eigvecs[:,res_1]),'b*')
    plt.plot(sp.sqrt(abs(2*problem.mass*Es[:,1])),abs(mb_eigvecs[:,res_1]),'g*')
    plt.figure(4)
    plt.clf()
    plt.plot(abs(mb_eigvecs[:,res_1]))
    plt.figure(5)
    plt.clf()
    plt.plot(abs(mb_eigvecs[dbl_res,:]),'ko')
    plt.plot(res_1, abs(mb_eigvecs[dbl_res,res_1]), 'go', markersize=8)
    plt.plot(res_2, abs(mb_eigvecs[dbl_res,res_2]), 'bo', markersize=8)

plot_shit(eigvals, mb_eigvals, mb_eigvecs)
plt.show()

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
