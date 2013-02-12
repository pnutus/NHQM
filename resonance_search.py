from imports import *
import basis_mom_space as mom
import problem_He5 as He5
import calculate_serial as calc
from wavefunction import gen_wavefunction
from normalize import normalize, norm

problem = He5.problem
steps = 100
lowest_energy = sp.empty(steps)
V0s = sp.linspace(-70, 0, steps)
l = 1
j = 1.5
order = 20

for (i, V0) in enumerate(V0s):
    problem.V0 = V0
    H = calc.hamiltonian(mom.H_element, args=(problem, l, j), order=order)
    energy, eigvecs = calc.energies(H)
    lowest_energy[i] = energy[0] * problem.eV_factor
    print V0
    
plt.plot(V0s, lowest_energy)
plt.show()
