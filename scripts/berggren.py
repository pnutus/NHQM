from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.bases import berggren_contour as berg
from nhqm.problems import He5
from nhqm.calculations import QM as calc

def gen_contour(tip_x, tip_y, k_max, points):
    """
    Generates a iosceles triangle contour 
    with tip at (tip_x, tip_y) which then 
    keeps going from (2*tip_x, 0) to (k_max, 0),
    """
    return berg.gen_berggren_contour(tip_x, - abs(tip_y), 
                k_max, points, first_zero_position = 2*tip_y)

problem = He5.problem   
order = 30
l = 1
j = 1.5
problem.V0 = -47.
step_size 

contour = gen_contour(0.17, 0.04, 2.5, order)
args = (problem, )
H = calc.hamiltonian(berg.H_complex_element, )