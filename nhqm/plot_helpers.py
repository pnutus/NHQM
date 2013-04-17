from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
from nhqm.QM_helpers import absq

def plot_contour(contour):
    points, _ = contour
    plt.plot(sp.real(points), sp.imag(points), 'o')

def plot_poles(energies, mass):
    ks = sp.sqrt(2*mass*energies)
    plt.plot(sp.real(ks), sp.imag(ks), 'o')
    
def plot_wavefunctions(contour, eigvecs):
    points, _ = contour
    plt.plot(sp.real(points), absq(eigvecs))
    
