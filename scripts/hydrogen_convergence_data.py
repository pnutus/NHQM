from __future__ import division
from imports import *
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
from nhqm.bases.gen_contour import gauss_contour
import plot_setup as plts
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os

def calc(ordermatrix, overwrite=False):
    
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_resho_path = "hydrogen_data/resho"
    rel_resm_path = "hydrogen_data/resm"
    abs_resho_file_path = os.path.join(script_dir, rel_resho_path)
    abs_resm_file_path = os.path.join(script_dir, rel_resm_path)
    
    
    
    try:
       with open(abs_resho_file_path) and open(abs_resm_file_path):
           if overwrite:
               raise IOError 
           resho = sp.loadtxt(abs_resho_file_path)
           resm = sp.loadtxt(abs_resm_file_path)
           return resho, resm
           
    except IOError:
        print "calculating"
    
        mom.integration_order = 20
        mom.integration_range = 10
        osc.integration_order = 60
        #ordermatrix=[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], [14,19,24,30,36,42,50]]

        ordermatrix_test=[[1,2], [14,19]]

        resho =[[],[]]
        resm=[[],[]]

        problem = H_atom
        omega =1
        problem.HO_omega = omega


        plts.plot_init(font_size=14,tick_size=11) #set font sizes.
        fig = plt.figure()

        for i, orderlist in enumerate(ordermatrix):
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
    
        sp.savetxt(abs_resho_file_path, resho)
        sp.savetxt(abs_resm_file_path, resm)        
            
     
        return resho, resm        

