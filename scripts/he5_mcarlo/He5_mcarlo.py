from __future__ import division
from nhqm.bases import harm_osc as osc, mom_space as mom, two_body_interaction as n_n
from nhqm.problems import H_atom, He5
from nhqm.QM_helpers import energies
import scipy as sp
from scipy.stats import uniform, norm
from nhqm.bases import gen_contour as cont
import matplotlib.pyplot as plt
from collections import namedtuple


def calc(vertices, order, num_points, reson = None, feed=None):    

    segs = len(vertices) - 1 
    QNums = namedtuple('qnums', 'l j J M E m')
    Q = QNums(l=1, j=1.5, J=0, M=0, 
              m=[-1.5, -0.5, 0.5, 1.5], 
              E=range(segs))
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
            print "example contour", contour
            print ""
            print ""
        if i % 500 == 0: 
            print i/ order, sum(res) / float(len(res))  
            
    #use previous results    
    if feed != None:    
        return (feed + sum(res) / float(len(res))  )/2
    else:
        return sum(res) / float(len(res))  
        
peak_x = 0.17
peak_y = -0.2*1j

He5.V0 = -70.

n_n.V0 = -1200
k_max = 5
delta = 0.05
vertices = [ 0, peak_x -delta + peak_y, peak_x +delta + peak_y, 2*peak_x, k_max]
vertices = [ 0, peak_x + peak_y, 2*peak_x, k_max]
num_points = 1 #per segment, only seems to work with 1, could it be that the result should be divided by num_points maybe? better, but not correct yet. it'd be such a shame if I had to think for myself
# ?
mom.integration_order = 20
mom.integration_range = 10
order = 10000
res = None
for i in xrange(1):
    res = calc(vertices,order, num_points, reson=(peak_x+peak_y), feed=res)
    print res, i    
    
 


"""conclusions, spillte carlo fuck up, doesn't converge and does so very slowly [sic]. if someone has the computer to do so, please let this run over night. yo"""  