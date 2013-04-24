from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour

import matplotlib.pyplot as plt

mom.integration_order = 20
mom.integration_range = 10
orderlist=[10,14,19,24,30,36,42,50]
#orderlist=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
resho =[]
resm=[]

problem = H_atom
omega =1
problem.HO_omega = omega

for order in orderlist:

    print ""
    Q = osc.QNums(l=0, j=.5, n=range(order))
    H = osc.hamiltonian(order, problem, Q)
    eigvals, eigvecs = energies(H)
    print osc.name, eigvals[0], problem.units
    resho.append(eigvals[0])
    
    contour = gauss_contour([0, 5], order)
    H = mom.hamiltonian(contour, problem, Q)
    eigvals, eigvecs = energies(H)
    print mom.name, eigvals[0], problem.units
    resm.append(eigvals[0])

    
plt.plot(orderlist,resho)
plt.plot(orderlist, resm)
plt.ylabel('some numbers')
plt.show() 