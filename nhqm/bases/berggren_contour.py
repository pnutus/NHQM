from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.special import sph_jn
import numpy as np
import mom_space as mom

name = "BerggrenContour"

@sp.vectorize    
def imag_integrand(r, k, k_prim, V, l, j):#takes complex p, p_prim and returns real Im of integrand
    j_l, _ = sph_jn(l, k*r)
    j_l_prim, _ = sph_jn(l, k_prim*r)
    return sp.imag( r**2 * j_l[-1] * j_l_prim[-1] * V(r, l, j) )
    

def H_complex_element(k, k_prim, problem, step_size, l = 0, j = .5):
    """step_size isn't entirelly correct but fr small peak amps should be close enough, 
    the contour generator should probably be able to handle this problem"""
    #print "hc " + str(p) + ", " + str(p_prim) 
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    
    real_integral, _ = fixed_quad(mom.real_integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j))

    imag_integral, _ = fixed_quad(imag_integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j)) 
                                                    
    return diagonal + 2 * k_prim**2 * step_size / sp.pi * (real_integral + 1j * imag_integral)
    
def real_berggren_wrapper(k_max, points):
    return gen_berggren_contour(0, 0, k_max, points, 0 )
    
def get_step_length( peak_re, peak_im, k_max, points, first_zero_position = -1 ):
    segment1 = sp.sqrt( peak_re ** 2 + peak_im ** 2)
    segment2 = sp.sqrt( peak_im ** 2 + (first_zero_position - peak_re) ** 2)
    
    length = segment1 + segment2 + (k_max - first_zero_position)
    
    return length / points    

def gen_berggren_contour(peak_re, peak_im, k_max, points, first_zero_position = -1 ):
    """generates a list with the coordinates for each point. """
    
    if first_zero_position == -1:
        first_zero_position = k_max
    
    segment1 = sp.sqrt( peak_re ** 2 + peak_im ** 2)
    segment2 = sp.sqrt( peak_im ** 2 + (first_zero_position - peak_re) ** 2)
    
    length = segment1 + segment2 + (k_max - first_zero_position)
    
    step_length = length / points
    
    #first segment, downward slope
    step_dw_re = peak_re / segment1 * step_length
    step_dw_im = peak_im / segment1 * step_length
    
    i=1
    while i * step_dw_im < peak_im:
        yield step_dw_re * i - 1j * step_dw_im * i
        i += 1
        
    #handling the peak
    peak_approx_re = step_dw_re * (i - 1)
    peak_approx_im = -step_dw_im * (i - 1)
    yield peak_approx_re + step_length + 1j * peak_approx_im
    
    #second segment, upward slope
    step_up_re = (first_zero_position - peak_re) / segment2 * step_length
    step_up_im = peak_im / segment2 * step_length
    
    """what if first_zero is k_max, might run out of points at the end of the upward slope, please implent a check for that."""
        
    j=1
    while peak_approx_im + step_up_im * j < 0:
        yield peak_approx_re + step_length + step_up_re * j + 1j * (peak_approx_im + step_up_im * j)
        j += 1
        
    #handling zero point
    zero_step_re = sp.sqrt( step_length ** 2 - (peak_approx_im + step_up_im * (j - 1)) ** 2 )

    zero_pos_approx = peak_approx_re + step_length + step_up_re * (j - 1) +1j * 0
    yield zero_pos_approx + zero_step_re + 1j * 0
    
    #we have now used (i - 1) + 1 + (j - 1) + 1 = i + j points; points -i -j remaining points.
    for k in xrange(points -i -j):
        yield zero_pos_approx + zero_step_re + step_length * k + 1j * 0
              
def gen_berggren_old_contour(peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1 ):              
    """if (peak_position == peak_neg_amplitude == first_zero_position == 0 or (peak_position == peak_neg_amplitude == 0 and first_zero_position ==  -1)):
        step0 = points
        step1 = step2 = 0
        
        for i in xrange(1, points+1):
            yield k_max / points * i
        
    else:"""
    step0 = int(np.ma.round( points * (peak_position / k_max) ) )
    step1 = int(np.ma.round( points * (first_zero_position - peak_position) / k_max ))
    step2 = int(np.ma.round( points * (k_max - first_zero_position) / k_max))
     
    segment_steps = [step0,step1,step2] 
     
    diff = points - np.sum(segment_steps)
    segment_steps[0] += diff
    
    re_step_size = [ 0.0, 0.0, 0.0]
    im_step_size = [ 0.0, 0.0, 0.0]
    
    #calculate step sizes in Re and Im for the downward slope

    if not peak_neg_amplitude == 0:
        re_step_size[0] = peak_position / segment_steps[0]
        im_step_size[0] = peak_neg_amplitude / segment_steps[0]
    
        re_step_size[1] = ( first_zero_position - peak_position) / segment_steps[1]
        im_step_size[1] = peak_neg_amplitude / segment_steps[1]
      
    
    if segment_steps[2] > 0:
        re_step_size[2] = (k_max - first_zero_position) / segment_steps[2]
    
    if not peak_neg_amplitude == 0:
        for point in xrange(1,segment_steps[0]+1):
            yield re_step_size[0]*point - 1j * im_step_size[0] * point
        
        for point in xrange(1,segment_steps[1]+1):
            yield peak_position + re_step_size[1] * point + \
            - 1j * peak_neg_amplitude + 1j * im_step_size[1] * point    

    for point in xrange(1, segment_steps[2]+1):
        yield first_zero_position + re_step_size[2] * point
              