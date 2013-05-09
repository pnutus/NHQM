from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import scipy as sp
import timeit

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

def absq(x):
    return x*sp.conjugate(x)

def logspace(min, max, order):
    return sp.around(sp.logspace(log(min), log(max), order))

problem = He5 
l = 1
j = 1.5
args = (problem, l, j)
k_max=5
order = 30

peak_x = 0.2
peak_y = 0.2
order_nr = 30
V0s = sp.linspace(-55, -40, order_nr)


#basis_function = mom.gen_basis_function(problem, l=l, j=j)
rmax=40
r = sp.linspace(1e-1, rmax, 2000)

plt.figure(2)
plt.title("Motion of the pole as V0 is varied from {3} MeV to {4} MeV \n  $l = {0},\, j = {1},\, N = {2}$".format(l, j, order, V0s[0], V0s[-1]), fontsize=20)
plt.xlabel('Re[k]', fontsize=20)
plt.ylabel('Im[k]', fontsize=20)
#print 'V0 =', V0,', peak x =', peak_x, 'peak y =', peak_y

k_res=sp.empty(order_nr, 'complex')
#eigvecs_all=sp.empty((order+2, order+2, cont_nrs), 'complex')
#eigvecs_res=sp.empty((order+2, order_nr), 'complex')
#absq_wavefunctions_res=sp.empty((2000, order_nr))
for m, V0 in enumerate(V0s):
    problem.V0=V0
    print order_nr-m, 'to go'
    contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
    ks, _ = contour
    Q = mom.QNums(l=1, j=1.5, k=range(len(contour[0])))
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = calc.energies(H)
    res = res_index(eigvecs)
    k_res[m]=sp.sqrt(2*problem.mass*eigvals[res])
    if sp.real(k_res[m])<10 ** -6:
        k_res[m] = abs(k_res[m]) * 1j
    #eigvecs_res[:,m] = eigvecs[:,res]
    #reswf = calc.gen_wavefunction(eigvecs[:,res], basis_function, contour)
    #reswf = calc.normalize(reswf, 0, 20, weight= lambda r: r**2)
    #absq_wavefunctions_res[:,m] = r**2 * absq(reswf(r))
print 'done'
plt.plot(sp.real(k_res), sp.imag(k_res), 'ko', markersize=2)

k=sp.sqrt(2*problem.mass*eigvals)
plt.plot(sp.real(k), sp.imag(k), 'ko', markersize=3)
plt.plot(sp.real(ks), sp.imag(ks), 'ro', markersize=3)

#plt.figure(1)
#plt.plot(r, absq_wavefunctions_res[:,0])


plt.show()
