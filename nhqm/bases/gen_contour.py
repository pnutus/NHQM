from __future__ import division
import scipy as sp
from scipy.special.orthogonal import p_roots

def gauss_contour(vertices, points_per_seg):
    """
    Generates a contour along the line segments between
    vertices using Gauss-Legendre quadrature.
    """
    
    num_segments = len(vertices) - 1
    points = weights = sp.empty(0, complex)
    for i in range(num_segments):
        (x, w) = p_roots(points_per_seg[i])
        a = vertices[i]
        b = vertices[i + 1]
        scaled_x = (x * (b - a) + (a + b))/2
        scaled_w = w * (b - a)/2
        points = sp.hstack((points, scaled_x))
        weights = sp.hstack((weights, scaled_w))
    return (points, weights)
    
def triangle_contour(peak_x, peak_y, k_max, order):
    return triangle_contour_explicit(peak_x, peak_y, k_max, 2*order, order)

def triangle_contour_explicit(peak_x, peak_y, k_max, trig_order, tail_order):
    vertices = [0, peak_x - 1j*peak_y, 2*peak_x, k_max]
    points_per_seg = [trig_order/2, trig_order/2, tail_order]
    return gauss_contour(vertices, points_per_seg)
    
# OLD CODE AHEAD!
    
@sp.vectorize
def triangle(x, peak_x, peak_y, k_max):
    if x < 2*peak_x:
        return peak_y * (abs(x/peak_x - 1) - 1)
    else:
        return 0

def naive_triangle_contour(peak_x, peak_y, k_max, points):
    xs = sp.linspace(0, k_max, points + 1)[1:]
    ys = triangle(xs, peak_x, peak_y, k_max)
    points = xs + 1j * ys
    weights = calculate_steps(points)
    return (points, weights)
    
def calculate_steps(contour):
    shifted = sp.roll(contour, 1)
    shifted[0] = 0
    return contour - shifted