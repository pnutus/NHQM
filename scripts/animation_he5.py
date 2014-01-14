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


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
ax.axis([-0.02, 0.7, -0.11, 0.01])


problem = He5 
order = 10*3
problem.V0 = -47.05
problem.Vso = -7.04
QNums = namedtuple('qnums', 'l j J M E m')
Q = QNums(l=1, j=1.5, J=0, M=0, 
          m=[-1.5, -0.5, 0.5, 1.5], 
          E=range(order))
k_max = 2.5

peak_x = 0.17


FPS = 24
film_length = 6
frames = film_length * FPS
peak_ys = sp.linspace(0, 0.1, frames)

def init():
    ks = ax.plot([],[])
    return ks, 

def animate(i):
    ax.cla()
    ax.axis([-0.1, 0.7, -0.11, 0.01])
    peak_y = peak_ys[i]
    contour = triangle_contour(peak_x, peak_y, k_max, order/3)
    points, weights = contour   
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    k = sp.sqrt(2*problem.mass*eigvals)
    mesh, = ax.plot(sp.real(points), sp.imag(points), 'og')
    ks, = ax.plot(sp.real(k), sp.imag(k), 'ob')
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    return ks, mesh, 

ani= anim.FuncAnimation(fig, animate, frames, init_func=init, interval=2, blit=False)
ani.save('res(contour).mp4', fps=FPS, writer=anim.FFMpegFileWriter())

#plt.show()
