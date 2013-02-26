from __future__ import division
import scipy as sp
from scipy.integrate import fixed_quad
from scipy.special import sph_jn
import numpy as np
import mom_space as mom

name = "Berggren"

@sp.vectorize
def real_integrand(*args):
    return sp.real( mom.integrand(*args) )

@sp.vectorize    
def imag_integrand(*args):
    return sp.imag( mom.integrand(*args) )
    

def H_element(k, k_prim, step, problem, l = 0, j = .5):
    diagonal = k**2 / (2 * problem.mass) * (k == k_prim)
    V = problem.potential
    
    real_integral, _ = fixed_quad(real_integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j))

    imag_integral, _ = fixed_quad(imag_integrand, 0, 10, \
                                n = 20, args=(k, k_prim, V, l, j)) 
                                                    
    return diagonal + 2 * k_prim**2 * step / sp.pi \
                        * (real_integral + 1j * imag_integral)

@sp.vectorize
def triangle(x, peak_x, peak_y, k_max):
    if x < 2*peak_x:
        return peak_y * (abs(x/peak_x - 1) - 1)
    else:
        return 0

def gen_simple_contour(peak_x, peak_y, k_max, points):
    xs = sp.linspace(0, k_max, points + 1)[1:]
    ys = triangle(xs, peak_x, peak_y, k_max)
    return xs + 1j * ys

def gen_contour(peak_re, peak_im, k_max, points, 
                        first_zero_position = None ):
    """generates a list with the coordinates for each point. """
    
    if first_zero_position == None:
        first_zero_position = k_max
    
    # calculate step length:
    
    segment1 = sp.sqrt( peak_re ** 2 + peak_im ** 2)
    segment2 = sp.sqrt( peak_im ** 2 + (first_zero_position - peak_re) ** 2)
    
    length = segment1 + segment2 + (k_max - first_zero_position)
    
    step_length = length / points
    
    # generate first segment, downward slope
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