from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour

import matplotlib.pyplot as plt


orderlist=[24,30,36,42,50,65,80,110,150]
res =[]

problem = H_atom
omega =4.7
problem.HO_omega = omega

for order in orderlist:

    print ""
    Q = osc.QNums(l=0, j=.5, n=range(order))
    H = osc.hamiltonian(order, problem, Q)
    eigvals, eigvecs = energies(H)
    print osc.name, eigvals[0], problem.units
    res.append(eigvals[0])

print len(orderlist),len(res)    
    
plt.plot(orderlist,res)
plt.ylabel('some numbers')
plt.show() 