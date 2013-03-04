from __future__ import division
from imports import *
import nhqm.calculations.QM as calc

def plot_contour(contour):
    points, _ = contour
    plt.plot(sp.real(points), sp.imag(points), 'o')

def plot_poles(energies, mass):
    ks = sp.sqrt(2*mass*energies)
    plt.plot(sp.real(ks), sp.imag(ks), 'o')
    
def plot_wavefunctions(contour, eigvecs):
    points, _ = contour
    plt.plot(sp.real(points), calc.absq(eigvecs))
    
