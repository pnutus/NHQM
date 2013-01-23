import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.interpolate import interpolate

def numerov(f, start, stop, dx, first, second):
    x = np.arange(start, stop, dx)
    psi = np.empty(len(x))
    dx2 = dx**2
    f1 = f(start)
    psi[0], psi[1] = first, second
    q0, q1 = psi[0]/(1 - dx2*f1/12), psi[1]/(1 - dx2*f1/12)
    for (i, ix) in enumerate(x):
        if i == 0: 
            continue
        f1 = f(ix)
        q2 = 2*q1 - q0 + dx2*f1*psi[i-1]
        q0, q1 = q1, q2
        psi[i] = q1/(1 - dx2*f1/12)
    return interpolate.interp1d(x, psi)
       
fun = numerov(lambda x: -1, 0, 10, 0.01, 0, 0.001)

x = np.arange(0,10,0.01)
y = fun(x)

plt.plot(x,y)
plt.show()