from __future__ import division
import scipy as sp
from scipy.special.orthogonal import p_roots

def gauss_contour(vertices, points_per_segment):
    """
    For each line between vertices
    """
    segment_count = len(vertices) - 1
    points = weights = sp.empty(0, complex)
    for i in range(segment_count):
        try:    (x, w) = p_roots(points_per_segment[i])
        except: (x, w) = p_roots(points_per_segment)
        a = vertices[i]
        b = vertices[i + 1]
        scaled_x = (x * (b - a) + (a + b))/2
        scaled_w = w * (b - a)/2
        points = sp.hstack((points, scaled_x))
        weights = sp.hstack((weights, scaled_w))
    return (points, weights)
        
def triangle_contour(peak_x, peak_y, k_max, points_per_segment):
    """
    Generates a triangular contour with a peak at (peak_x, -peak_y)
    ending at k_max, with a certain number of points per segment.
    
          _ _ _ _ k_max
    \    /
     \  /
      \/
     peak
    
    """
    return triangle_contour_explicit(peak_x, peak_y, k_max, 
                                    2*points_per_segment, points_per_segment)

def triangle_contour_explicit(peak_x, peak_y, k_max, 
                              triangle_point_count, tail_point_count):
    """
    Generates a triangular contour with 
    """
    vertices = [0, peak_x - 1j*peak_y, 2*peak_x, k_max]
    points_per_segment = [ triangle_point_count/2
                         , triangle_point_count/2
                         , tail_point_count
                         ]
    return gauss_contour(vertices, points_per_segment)
    
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