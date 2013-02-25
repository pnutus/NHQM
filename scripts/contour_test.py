from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom, harm_osc as osc, berggren as berg
from nhqm.problems import He5, H_atom
from nhqm.calculations import QM as calc
import time
import numpy as np
import re

def find(v,x):
    for i in xrange(len(v)):
        if (v[i] == x):
            print i
            return i
    print str(-1)        
    return -1
                 
def kontroll_samma_kontur(order, k_max, peak_pos, peak_amp, zero_point, l=0, j=0.5):
    #order = 10
    problem = He5.problem
    basis = mom
    #k_max = 3   #only k_max = order works with the current contour implementation...
    step_size = k_max / order #only step_size = 1 works with the current contour implementation...
    #l = 0
    #j = 0.5
    r = sp.linspace(0, 10, 200)

    kontur_a = berg.gen_contour( peak_pos, peak_amp, k_max, order, zero_point )
    kontur_b = berg.gen_contour( peak_pos, peak_amp, k_max, order, zero_point )

    a = 0
    kont_hamilton = sp.empty( (order,order), complex ) # complex elements?
    for k in kontur_a:
        b = 0
        for k_prim in kontur_b:
            kont_hamilton[a,b] = berg.H_element(k, k_prim, He5.problem, step_size, l, j)
            b += 1
        a += 1
        kontur_b = berg.gen_contour( peak_pos, peak_amp, k_max, order, zero_point )  
    
    
    H = calc.hamiltonian(basis.H_element, args=(problem, step_size, l, j), order=order)

    return kont_hamilton , H

def plott_kontur( peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1  ):
        return plot_contour( peak_position, peak_neg_amplitude, k_max, points, first_zero_position )
    
def plot_contour( peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1 ):    
    cont = berg.gen_contour( peak_position, peak_neg_amplitude, k_max, points, first_zero_position )
    res = []
    for i in cont:
         res.append(i)
         
    re = sp.real(res)
    im = sp.imag(res)
    
    plt.xlabel('Re')
    plt.ylabel('Im')
    #plt.title('Berggren kontur')
    

    plt.plot(re, im)

    lin = sp.linspace(0,k_max)
    noll = sp.zeros((len(lin),1))
    
    plt.plot(lin, noll, ls ='dashed')
    plt.show()
    
    return res
    #print res
    
def plot_eigenvectors(order, k_max, peak_pos, peak_amp, zero_point, l=0, j=0.5):
    #this should only take amatrix and be a general method, but fuck generality
    h_c, _ = kontroll_samma_kontur(order, k_max, peak_pos, peak_amp, zero_pos, l, j)
    #print hc
    #print h
    
    plt.subplot(2,1,1)
    fig1 = plt.gcf()
    plt.xlabel('Contour point index')
    plt.ylabel('abs(phi)')
    plotttitle= "Eigenvectors (abs(phi)) for He5 problem (l = " + str(l) +" j = " + str(j) + ") in berggren basis \n with a contour triangle  0, (" + str(peak_pos) + ", -" + str(peak_amp) + "j), (" + str(zero_pos) + ",0j) and then a line to k_max = " + str(k_max) + "; to order " + str(order)
    plt.title(plotttitle)

    eigval_c, eigvec_c = sp.linalg.eig(h_c)
    
    vec_max = sp.empty((order,1), complex)
    x = sp.linspace(0,order,order)
    for i in xrange(0,order):
        temp=np.abs(eigvec_c[:,i])
        plt.plot(x,temp)
        vec_max[i] = max(temp)
    min_index = find( vec_max, min(vec_max) )
    
    plt.plot(x,np.abs(eigvec_c[:,min_index]), lw = 3)
        
    plt.subplot(2,1,2)
    #peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1 
    plot_contour(peak_pos, peak_amp, k_max, order, zero_pos)
    
    filetitle = "He5 problem (l = " + str(l) +" j = " + str(j) + ") 0, (" + str(peak_pos) + ", -" + str(peak_amp) + "j), (" + str(zero_pos) + ",0j) k_max = " + str(k_max) + "; to order " + str(order)
    filetitle = re.sub('[.]', ',', filetitle)
    
    
    #figpath = "/Users/Spill/Documents/kandidatarbete/konturplotar/" + filetitle
    #fig1.savefig(figpath)
    plt.show()    
    
    return fig1
if __name__ == '__main__':
    
    #peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1 
    
    #plott_kontur( 1, .1, 10, 100, 5)
    
    order = 100
    k_max = 5
    peak_pos = 1
    peak_amp = 1.5
    zero_pos = 3
    l=0
    j=0.5
    #print plot_contour( peak_pos, peak_amp, k_max, order, zero_pos )
    plot_eigenvectors( order, k_max, peak_pos, peak_amp, zero_pos, l=0, j=0.5 )