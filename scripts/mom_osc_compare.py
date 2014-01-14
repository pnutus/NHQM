from __future__ import division
from imports import *
from nhqm.bases.mom import mom_space as mom, harm_osc as osc
from nhqm.problems import He5, H_atom
from nhqm.calculations import QM as calc

problems = [He5, H_atom]
bases = [mom, osc]
order = 30
k_max = 5
step_size = k_max / order
l = 0
j = 0.5
r = sp.linspace(0, 10, 200)
    
plt.suptitle("Radial probability distribution "
          "with $l = {1}$, $j = {2}$ "
          "and a ${3} \\times\\,{3}$ matrix"
          .format("", l, j, order))

for (i, problem) in enumerate(problems):
    print problem.name
    plt.subplot(1, 2, i+1)
    plt.title(problem.name)
    plt.ylabel(r'$r^2|\Psi|^2$')
    plt.xlabel(r'$r$')
    plt.show(block=False)
    
    
    for (i, basis) in enumerate(bases):
        if basis is osc:
            omega = osc.optimal_osc_freq(problem, l, j)
            H = calc.hamiltonian(osc.H_element, args=(problem, omega, l, j), order=order)
            basis_function = basis.gen_basis_function(problem, omega, l=l, j=j)
        else:
            H = calc.hamiltonian(basis.H_element, args=(step_size, problem, l, j), order=order)
            basis_function = basis.gen_basis_function(step_size, problem, l=l, j=j)
        energy, eigvecs = calc.energies(H)
        lowest_energy = energy[0] * problem.eV_factor
        print basis.name, "lowest energy:", lowest_energy, "eV"
        wavefunction = calc.gen_wavefunction(eigvecs[:, 0], basis_function)
        wavefunction = calc.normalize(wavefunction, 0, 10, weight=lambda r: r**2)
        
        plt.plot(r, r**2 * absq(wavefunction(r)), label = basis.name)
        plt.legend()
        plt.draw()

plt.show() # pause when done