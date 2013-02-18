from imports import *
import basis_mom_space as mom
import problem_He5 as He5
import calculate_serial as calc
from wavefunction import gen_wavefunction
from normalize import normalize, norm

problem = He5.problem
problem.V0 = -52.3
order = 20
l = 1
j = 1.5

def absq(x):
    return x*sp.conjugate(x)

H = calc.hamiltonian(mom.H_element, args=(problem, l, j), order=order)
energy, eigvecs = calc.energies(H)
basis_function = mom.gen_basis_function(problem, l=l, j=j)
wavefunction = gen_wavefunction(eigvecs[:,0], basis_function)
wavefunction = normalize(wavefunction, 0, 10, weight= lambda r: r**2)


r = sp.linspace(1e-1, 10, 200)

plt.title("Potential, ground energy (and wavefunction, not to scale) \n $l = {0}, j = {1}, V_0 = {2}$MeV")
plt.plot(r, l*(l + 1)/ r**2 + problem.potential(r, l, j))
plt.plot(r, 50*r**2 * absq(wavefunction(r)) + energy[0])
plt.plot(r, sp.zeros(len(r)) + energy[0])
plt.axes([0, 10, -0.05, 0.05])
plt.ylabel(r'$r^2|\Psi(r)|^2$')
plt.xlabel(r'$r$ / fm')
plt.show()