from imports import *
from scipy.interpolate import interpolate

def numerov(f, start, stop, dx, first, second):
    x = np.linspace(start, stop, abs(start - stop)/dx )
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
       
# def V(r):
#     return - 1. / r
#        
# sol = numerov(lambda r: -1, 1e-10, 100, 0.01, 0, 1)
# 
# x = np.arange(1e-10,100,0.01)
# y = sol(x)
# 
# plt.plot(x,y)
# plt.show()