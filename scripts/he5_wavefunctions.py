from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import triangle_contour, naive_triangle_contour, gauss_contour
from nhqm.plot_helpers import *
from numpy.linalg import norm
from collections import namedtuple

#Plots contour, mom space wavefunctions (normed |phi(|k|)/ksqrt(w)|), coordinate space wavefunctions


def res_index():
    result = 0
    global_min = 1
    for i in range(order):
        temp_max = max( abs(eigvecs[:,i]/sp.sqrt(weights))/norm(abs(eigvecs[:,i]/sp.sqrt(weights))) )
        if temp_max < global_min:
            global_min = temp_max
            result = i
    return result

def prob_wf(index):
    wf = mom.gen_wavefunction(eigvecs[:,index], Q, contour)
    return r_order / rmax * r**2 * absq(wf(r))/ norm(r*wf(r))**2

def plot_mom_wf(index, c = '', lw = 1):
    wf = abs(eigvecs[:,index]/(sp.sqrt(weights)))/norm(abs(eigvecs[:,index]/(sp.sqrt(weights))))
    plt.plot(sp.real(points), wf ** 2, c, linewidth = lw)

problem = He5 
order = 15*3
problem.V0 = -47.05
problem.Vso = -7.04
peak_x = 0.35
peak_y = 0.4
k_max = 2.5



contour = triangle_contour(peak_x, peak_y, k_max, order/3)
#contour = gauss_contour([0, k_max], order)
#contour = naive_triangle_contour(peak_x, peak_y, k_max, order)
points, weights = contour

QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(order))



H = mom.hamiltonian(contour, problem, Q)
eigvals, eigvecs = energies(H)


rmax = 20
r_order = 300
r = sp.linspace(0, rmax, r_order)

plt.figure(1)
plt.clf()
plt.title('$k=sqrt(2mE)$ and mesh points \n $V0 = {0}$, $N = {1}$'.format(problem.V0, order))
plot_poles(eigvals, problem.mass)
plt.plot(sp.real(points), sp.imag(points), 'og')
plt.axis([0,2*peak_x, -1.1*peak_y, 0])

plt.figure(2)
plt.clf()
plt.title('abs($\phi(k)$/(k sqrt(weights)) \n $V0 = {0}$, $N = {1}$'.format(problem.V0, order))
for m in range(order):
    plot_mom_wf(m)

res = res_index()
plot_mom_wf(res, 'k', 3)
plot_mom_wf(res-1, 'g', 3)

plt.figure(3)
plt.clf()
plt.title('$r^2|\psi|^2$ \n $V0 = {0}$, $N = {1}$'.format(problem.V0, order))
plt.plot(r, prob_wf(res_index()),'k', linewidth=3)
plt.plot(r, prob_wf(res_index()-1),'g', linewidth=3)

plt.figure(4)
plt.clf()
plt.title('$|\phi(|k|)$ along the contour \n $V0 = {0}$, $N = {1}$'.format(problem.V0, order))
plot_mom_wf(res, 'k', 4)

plt.show()

