def custom_integrate(function, contour, args=()):
    sum = 0
    for point, weight in zip(*contour):
        sum += weight * function(point, *args)
    return sum
    
def f(x,y,z):
    return x*x