from imports import *
import basis_mom_space
import basis_harm_osc
import problem_He5
import calculate_serial as calc

problem = problem_He5.problem

l = 0
j = .5

H_mom = calc.hamiltonian(basis_mom_space.H_element, args=(problem, l, j))
print "MomSpace lowest energy:", calc.energies(H_mom)[0] * problem.eV_factor, "eV"

omega = basis_harm_osc.optimal_osc_freq(problem, l, j)
H_osc = calc.hamiltonian(basis_harm_osc.H_element, args=(problem, omega, l, j))
print "HarmOsc lowest energy:", calc.energies(H_osc)[0] * problem.eV_factor, "eV"