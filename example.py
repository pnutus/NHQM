from imports import *
import basis_mom_space as mom
import basis_harm_osc as osc
import problem_He5 as He5
import problem_H_atom as H_atom
import calculate_serial as calc
from wavefunction import gen_wavefunction
from scipy import linalg


problems = [H_atom.problem, He5.problem]
order = 20
l = 0
j = .5

def absq(x):
    return x*sp.conjugate(x)

for (i, problem) in enumerate(problems):
    print problem.name
    
    # Lowest energy using momentum basis
    H_mom = calc.hamiltonian(mom.H_element, args=(problem, l, j), order=order)
    lowest_energy = calc.energies(H_mom)[0] * problem.eV_factor
    print "MomSpace lowest energy:", lowest_energy, "eV"
    
    # Lowest energy using harmonic oscillator basis
    omega = osc.optimal_osc_freq(problem, l, j)
    H_osc = calc.hamiltonian(osc.H_element, args=(problem, omega, l, j), order=order)
    lowest_energy = calc.energies(H_osc, hermitian=True)[0] * problem.eV_factor
    print "HarmOsc lowest energy:", lowest_energy, "eV"
    
    _, eigvecs = linalg.eigh(H_osc)
    basis_function = osc.gen_basis_function(problem, omega, l, j)
    wavefunction = gen_wavefunction(eigvecs[:,0], basis_function)
    r = sp.linspace(0, 10, 100)
    plt.figure(i)
    plt.title("{0} radial wave function\n"
              "with $l = {1}$, $j = {2}$ "
              "and a ${3} \\times {3}$ matrix"
              .format(problem.name, l, j, order))
    plt.ylabel('$|\Psi|^2$')
    plt.xlabel('$r$')
    plt.plot(r, r**2 * absq(wavefunction(r)))
    plt.show(block=False)
    
plt.show() # pause when done