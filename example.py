from imports import *
import basis_mom_space as mom
import basis_harm_osc as osc
import problem_He5 as He5
import problem_H_atom as H_atom
import calculate_serial as calc
from wavefunction import gen_wavefunction
from normalize import normalize, norm

problems = [He5.problem, H_atom.problem,]
bases = [mom, osc]
order = 20
l = 1
j = 0.5
r = sp.linspace(0, 10, 200)

def absq(x):
    return x*sp.conjugate(x)
    
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
            H = calc.hamiltonian(basis.H_element, args=(problem, l, j), order=order)
            basis_function = basis.gen_basis_function(problem, l=l, j=j)
        energy, eigvecs = calc.energies(H)
        lowest_energy = energy[0] * problem.eV_factor
        print basis.name, "lowest energy:", lowest_energy, "eV"
        wavefunction = gen_wavefunction(eigvecs[:, 0], basis_function)
        wavefunction = normalize(wavefunction, 0, 10, weight=lambda r: r**2)
        
        plt.plot(r, r**2 * absq(wavefunction(r)), label = basis.name)
        plt.legend()
        plt.draw()

plt.show() # pause when done