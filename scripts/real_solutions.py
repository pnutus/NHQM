from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.bases import many_body as mb
from nhqm.problems import He5
from nhqm.bases.gen_contour import triangle_contour, naive_triangle_contour, gauss_contour
from nhqm.QM_helpers import gen_wavefunction, energies
from numpy.linalg import norm

def absq(x):
    return x*sp.conjugate(x)


problem = He5
k_max = 3
order = 100

peak_x = 0
peak_y = 0
contour = gauss_contour((0, k_max), order)
#contour = naive_triangle_contour(peak_x, peak_y, k_max, order)
ks, _ = contour
Q = mb.QNums(l=1, J=0, M=0, j=1.5, 
             m=[-1.5, -0.5, 0.5, 1.5], 
             k=ks)

problem.V0 = -52
H = mom.hamiltonian(contour, problem, Q)
eigvals_1, eigvecs_1 = energies(H)

problem.V0=-47
H = mom.hamiltonian(contour, problem, Q)
eigvals_2, eigvecs_2 = energies(H)


basis_function = mom.gen_basis_function(problem, Q.l, Q.j)
rmax = 100
r_order = 500
r = sp.linspace(1e-1, rmax, r_order)

def sqrd_wf(eigvec):
    wf = gen_wavefunction(eigvec, basis_function, contour=contour)
    wf=wf(r)
    return r_order / rmax * r**2 * absq(wf) / norm(r*wf)**2 

def plotwf(eigvec, energy, color='k'):
    wf = sqrd_wf(eigvec)
    plt.plot((0, rmax), (energy, energy),'r',linewidth=2)
    plt.plot(r, wf + energy, color, linewidth=3)
    return _



plt.figure(1)
plt.clf()
plt.title("Eigensolution wavefunctions. \n  $V0 = {0},\, N = {1}$".format(problem.V0, order), fontsize=20)
plt.xlabel('r', fontsize=20)
plt.ylabel('$r^2|R(r)|^2$', fontsize=20)

problem.V0=-52
V = Q.l * (Q.l + 1) / r**2 + problem.potential(r, Q.l, Q.j)
plt.plot(r, Q.l * (Q.l + 1) / r**2 + V, '--k', linewidth=3)

plotwf(eigvecs_1[:,10], eigvals_1[10])
plotwf(eigvecs_1[:,14], eigvals_1[14])
plotwf(eigvecs_1[:,18], eigvals_1[18])
plotwf(eigvecs_1[:,7], eigvals_1[7])

plt.figure(2)
plt.clf()
plt.title("Eigensolution wavefunctions. \n  $V0 = {0},\, N = {1}$".format(problem.V0, order), fontsize=20)
plt.xlabel('r', fontsize=20)
plt.ylabel('Energy [MeV] / $r^2|R(r)|^2$', fontsize=20)

problem.V0=-47
V = Q.l * (Q.l + 1) / r**2 + problem.potential(r, Q.l, Q.j)
plt.plot(r, Q.l * (Q.l + 1) / r**2 + V, '--k', linewidth=3)

plotwf(eigvecs_2[:,10], eigvals_2[10])
plotwf(eigvecs_2[:,14], eigvals_2[14])
plotwf(eigvecs_2[:,18], eigvals_2[18])
plotwf(eigvecs_2[:,7], eigvals_2[7])


plt.figure(3)
plt.clf()
plt.plot(abs(eigvecs_1))
plt.plot(abs(eigvecs_1[:,7]),'k', linewidth=3)

plt.figure(4)
plt.clf()
plt.plot(abs(eigvecs_2))
plt.plot(abs(eigvecs_2[:,7]),'k', linewidth=3)


plt.show()

