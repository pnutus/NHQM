from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.bases import many_body as mb
from nhqm.problems import He5
from nhqm.bases.gen_contour import triangle_contour, naive_triangle_contour, gauss_contour
from nhqm.QM_helpers import gen_wavefunction, energies, absq, normalize
from numpy.linalg import norm


problem = He5
k_max = 3
order = 100

peak_x = 0
peak_y = 0
contour = gauss_contour((0, k_max), order)
#contour = naive_triangle_contour(peak_x, peak_y, k_max, order)
ks, _ = contour
Q = mom.QNums(l=1, j=1.5, k=ks)

problem.V0 = -52
H = mom.hamiltonian(contour, problem, Q)
eigvals_1, eigvecs_1 = energies(H)

problem.V0=-50
H = mom.hamiltonian(contour, problem, Q)
eigvals_2, eigvecs_2 = energies(H)


rmax = 100
r_order = 500
r = sp.linspace(1e-1, rmax, r_order)

def sqrd_wf(eigvec):
    wf = mom.gen_wavefunction(eigvec, contour, Q)
    return r_order / rmax * r**2 * absq(wf(r))*10/ norm(r*wf(r))**2 

def plotwf(eigvec, energy, color='k'):
    wf = sqrd_wf(eigvec)
    plt.plot((0, rmax), (energy, energy),'r',linewidth=2)
    plt.plot(r, wf + energy, color, linewidth=3)



#plt.figure(1)
#plt.clf()
#plt.title("Eigensolution wavefunctions. \n  $V0 = {0},\, N = {1}$".format(problem.V0, order), fontsize=20)
#plt.xlabel('r', fontsize=20)
#plt.ylabel('$r^2|R(r)|^2$', fontsize=20)

#plotwf(eigvecs_1[:,8], eigvals_1[8])
#plotwf(eigvecs_1[:,17], eigvals_1[17])
#plotwf(eigvecs_1[:,19], eigvals_1[19])
#plotwf(eigvecs_1[:,20], eigvals_1[20])

#plt.figure(2)
#plt.clf()
#plt.title("Eigensolution wavefunctions. \n  $V0 = {0},\, N = {1}$".format(problem.V0, order), fontsize=20)
#plt.xlabel('r', fontsize=20)
#plt.ylabel('Energy [MeV] / $r^2|R(r)|^2$', fontsize=20)

#plotwf(eigvecs_2[:,13], eigvals_2[13])
#plotwf(eigvecs_2[:,17], eigvals_2[17])
#plotwf(eigvecs_2[:,19], eigvals_2[19])
#plotwf(eigvecs_2[:,20], eigvals_2[20])

 
#plt.figure(3)
#plt.clf()
#plt.plot(abs(eigvecs_1))
#plt.plot(abs(eigvecs_1[:,8]),'k', linewidth=3)
 
#plt.figure(4)
#plt.clf()
#plt.plot(abs(eigvecs_2))
#plt.plot(abs(eigvecs_2[:,13]),'k', linewidth=3)


wf1_res = sqrd_wf(eigvecs_1[:,8])
wf1_1 = sqrd_wf(eigvecs_1[:,17])

wf2_res = sqrd_wf(eigvecs_2[:,13])
wf2_1 = sqrd_wf(eigvecs_2[:,17])


fig1=plt.figure(5)
fig1.clf()
ax1=fig1.add_subplot(121)
ax2=fig1.add_subplot(122)

ax1.plot(r, wf1_res, 'k', linewidth=3)
ax1.plot(r, wf1_1, 'k', linewidth=3) 

ax2.plot(r, wf_res, 'k', linewidth=3)
ax2.plot(r, wf_1, 'k', linewidth=3) 



ax1.axis([0, 100, 0, 1.6])
ax1.get_yaxis().set_ticks([])
ax1.set_xlabel('$r$ [fm]', fontsize=24)
ax1.set_ylabel('$r^2|R(r)|^2$', fontsize=24)
ax1.set_title('$V0 = -52$', fontsize=24)
plt.setp(ax1.get_xticklabels(), fontsize=24)

ax2.axis([0, 100, 0, 1.6])
ax2.get_yaxis().set_ticks([])
ax2.set_xlabel('$r$ [fm]', fontsize=24)
ax2.set_ylabel('$r^2|R(r)|^2$', fontsize=24)
ax2.set_title('$V0 = -47$', fontsize=24)
plt.setp(ax2.get_xticklabels(), fontsize=24)

plt.show()

