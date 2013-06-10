from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.QM_helpers import energies
from nhqm.bases.contours import triangle_contour, naive_triangle_contour, gauss_contour
from nhqm.plot_helpers import *
from numpy.linalg import norm
from collections import namedtuple
import matplotlib.animation as anim


fig = plt.figure(2)
fig.clf()
ax = fig.add_subplot(111)
ax.axis([-0.02, 0.7, -0.15, 0.01])


problem = He5 
order = 20*3
V0max = 60
V0min = 35
problem.Vso = -7.04
QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(order))
k_max = 2.5

peak_x = 0.25
peak_y = 0.15

contour = triangle_contour(peak_x, peak_y, k_max, order/3)
points, weights = contour   


FPS = 40
film_length = 20
frames = film_length * FPS
V0s = sp.linspace(-V0max, -V0min, frames)

def init():
    ks = ax.plot([],[])
    return ks, 

def animate(i):
    ax.cla()
    ax.axis([-0.01, 1, -0.2, 0.25])
    ax.set_xlabel('Re k / [fm$^{-1}$]', fontsize=20)
    ax.set_ylabel('Im k / [fm$^{-1}$]', fontsize=20)
    problem.V0 = V0s[i]
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    
    k = sp.sqrt(2*problem.mass*eigvals)
    for a, b in enumerate(k):
        if sp.real(b)<1e-3:
            k[a]=sp.real(b)+abs(sp.imag(b))*1j
    mesh, = ax.plot(sp.real(points), sp.imag(points), 'o',color='0.5')
    ks, = ax.plot(sp.real(k), sp.imag(k), 'or')
    text = ax.text(0.5, 0.17, 'V0 = %.0f MeV'%(problem.V0), fontsize=24)
    return ks, mesh, text

ani= anim.FuncAnimation(fig, animate, frames, init_func=init, interval=2, blit=False)
ani.save('polsvep_anim_HD.mp4', fps=FPS)

plt.show()
