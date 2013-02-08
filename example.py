from imports import *
import basis_mom_space as mom
import basis_harm_osc as osc
import problem_He5 as He5
import problem_H_atom as H_atom
import calculate_serial as calc


problems = [H_atom.problem, He5.problem]
order = 20
l = 0
j = .5

for problem in problems:
    print problem.name
    
    H_mom = calc.hamiltonian(mom.H_element, args=(problem, l, j), order=order)
    lowest_energy = calc.energies(H_mom)[0] * problem.eV_factor
    print "MomSpace lowest energy:", lowest_energy, "eV"

    omega = osc.optimal_osc_freq(problem, l, j)
    H_osc = calc.hamiltonian(osc.H_element, args=(problem, omega, l, j), order=order)
    lowest_energy = calc.energies(H_osc, hermitian=True)[0] * problem.eV_factor
    print "HarmOsc lowest energy:", lowest_energy, "eV"