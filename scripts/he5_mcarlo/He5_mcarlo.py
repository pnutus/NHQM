from __future__ import division
from nhqm.bases import harm_osc as osc, mom_space as mom
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
import scipy as sp
from scipy.stats import uniform, norm
from nhqm.bases import gen_contour as cont
import matplotlib.pyplot as plt

def calc(vertices, order, num_points, reson = None, feed=None):    
    mom.integration_order = 20
    mom.integration_range = 10

    segs = len(vertices) - 1
    Q = osc.QNums(l=0, j=.5, n=range(segs))
    randoms = uniform.rvs(size=(order, segs, num_points))
        
    res=[]

    problem = He5

    for i in xrange(order):
        contour = []
        points = []
        for j in xrange(segs):
            randoms[i][j].sort()
            for k in xrange(num_points): #num_points: probably a bad idea
                point = vertices[j] + (vertices[j+1] - vertices[j])*randoms[i][j][k]
                points.append(point)
                if j == 1 and k == 0:
                    points.append(reson)
        
        weights = cont.calculate_steps(points)
        contour = points, weights
        H = mom.hamiltonian(contour, problem, Q)
        eigvals, eigvecs = energies(H)
        res.append(eigvals[0])
        
        if i == 0:
            print contour
        if i % 500 == 0: 
            print i/ order, sum(res) / float(len(res))  
            
    #use previous results    
    if feed != None:    
        return (feed + sum(res) / float(len(res))  )/2
    else:
        return sum(res) / float(len(res))  
        
peak_x = 0.17
peak_y = -0.07*1j
k_max = 2.5
delta = 0.05
vertices = [ 0, peak_x -delta + peak_y, peak_x +delta + peak_y, 2*peak_x, k_max]
vertices = [ 0, peak_x + peak_y, 2*peak_x, k_max]
num_points = 2 #per segment, only seems to work with 1, could it be that the result should be divided by num_points maybe? better, but not correct yet. it'd be such a shame if I had to think for myself
# ?
mom.integration_order = 20
mom.integration_range = 10
order = 10000
res = None
for i in xrange(1):
    res = calc(vertices,order, num_points, reson=(peak_x+peak_y), feed=res)
    print res, i    
    
"""calculating
[[(-45.827216921250169+0.0031556789657023951j), (-35.19737169386039+0.00039054401046521446j), (-23.982400963077929-1.9178313356026469e-06j), (-24.826392768520193-4.4877508943058776e-10j), (-24.937558214753313+1.1102491368587253e-11j), (-24.913746184094386-2.6477506591142258e-14j), (-24.913117449957419-1.3906593612060391e-15j), (-24.913275768411108+1.2682017131532579e-15j), (-24.913286332740409-1.6601775356638704e-15j), (-24.913285224535802-8.1357415683585126e-16j)], [(-24.913285224535802-8.1357415683585126e-16j), (-24.913285072081791+6.6626824050227484e-17j), (-24.913285072081308-2.6594775427056549e-16j), (-24.913285072081319-1.7091075140844625e-16j), (-24.913285072081408+2.5602631676850566e-15j), (-24.913285072081369-1.2055249999533334e-15j),
^^^He5_conv_plot calculations for the triangle contour for different (rising) orders

THIS IS WHAT WE ARE HOPING FOR (-24.913285072081443+2.4589456699190438e-15j)"""  


"""conclusions, spillte carlo fuck up, doesn't converge and does so very slowly [sic]. if someone has the computer to do so, please let this run over night. yo"""  