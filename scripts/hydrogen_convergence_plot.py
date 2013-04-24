from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour
import plot_setup as plts
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

mom.integration_order = 20
mom.integration_range = 10
osc.integration_order = 60
orderpappa=[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], [14,19,24,30,36,42,50]]
resho =[[],[]]
resm=[[],[]]

problem = H_atom
omega =1
problem.HO_omega = omega


plts.plot_init(font_size=14,tick_size=11) #set font sizes.
fig = plt.figure()

#fixes subplots, and size
gs = gridspec.GridSpec(1,2,
                       width_ratios=[1,2]
                       )
ax2 = plt.subplot(gs[0])
ax3 = plt.subplot(gs[1])
ax_list = [ax2, ax3]


ax2.set_ylabel('energy [Hartree]')
plot_title=plt.title('Convergence of the Hydrogen groundstate for HO, mom')
#ad hoc solutions to title position
plot_title.set_y(1.03)
plot_title.set_x(0.13)


for i, orderlist in enumerate(orderpappa):
    #print orderlist
    for order in orderlist:

        print ""
        Q = osc.QNums(l=0, j=.5, n=range(order))
        H = osc.hamiltonian(order, problem, Q)
        eigvals, eigvecs = energies(H)
        #print osc.name, eigvals[0], problem.units
        resho[i].append(eigvals[0])
    
        contour = gauss_contour([0, 5], order)
        H = mom.hamiltonian(contour, problem, Q)
        eigvals, eigvecs = energies(H)
        #print mom.name, eigvals[0], problem.units
        resm[i].append(eigvals[0])
    
    ax = ax_list[i]
    l1, l2 = ax.plot(orderlist, resho[i], 'bo-', orderlist, resm[i], 'go-')
    ax.set_xlabel('order')

    if i==1:
        ax.set_ylim([-0.505, -0.485])
        ax.legend( (l1, l2),
                ('Harmonic Oscillator', 'Momentun Space'),
                'lower right')
plt.show() 