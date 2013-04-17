def gauss_contour(vertices, order):
    """
    Generates a contour along the line segments between
    vertices using Gauss-Legendre quadrature.
    """
    (x, w) = p_roots(order)
    num_segments = len(vertices) - 1
    points = weights = sp.empty(0, complex)
    for i in range(num_segments):
        a = vertices[i]
        b = vertices[i + 1]
        scaled_x = (x * (b - a) + (a + b))/2
        scaled_w = w * (b - a)/2
        points = sp.hstack((points, scaled_x))
        weights = sp.hstack((weights, scaled_w))
    return (points, weights)
    
def triangle_contour(peak_x, peak_y, k_max, order):
    vertices = [0, peak_x - 1j*peak_y, 2*peak_x, k_max]
    return gauss_contour(vertices, order)
    
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