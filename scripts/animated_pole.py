from __future__ import division
from imports import *
from nhqm.plot_helpers import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from helpers import progress
from nhqm.bases.gen_contour import triangle_contour

problem = He5  
order = 30
peak_x = 0.17
peak_y = 0.07
k_max = 2.5
steps = 5

Q = mom.QNums(l=1, j=1.5, k=range(order))

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

contour = triangle_contour(peak_x, peak_y, k_max, order/3)
V0s = sp.linspace(-53, -51, steps)
poles = sp.empty(steps, complex)

for (i, V0) in enumerate(progress(V0s)):
    problem.V0 = V0
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    idx = res_index(eigvecs)
    res_E = eigvals[idx]
    pole = sp.sqrt(2*problem.mass*res_E)
    if sp.real(pole) < 1e-4 and abs(sp.imag(pole)) > 0.02:
        pole = 1j*abs(pole)
    poles[i] = pole

plot_contour(contour)
# plot_poles(eigvals, problem.mass)
plt.plot(sp.real(poles), sp.imag(poles), 'ro')
plt.show()