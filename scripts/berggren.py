from __future__ import division
from imports import *
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import timeit

def res_index(eigvecs):
    maxes = map(max, abs(eigvecs.T))
    return maxes.index(min(maxes))

def absq(x):
    return x*sp.conjugate(x)

problem = He5.problem   
l = 1
j = 1.5
args = (problem, l, j)
k_max=2.5
order = 9

cont_nrs = 3
peak_x = 0.173
peaks_y = sp.linspace(.05, .1, cont_nrs)
V0 = -47.
problem.V0 = V0


rmax=50
r = sp.linspace(1e-1, rmax, 2000)

basis_function = mom.gen_basis_function(problem, l=l, j=j)






plt.figure(1)
#plt.title("The figure shows the momenta k=sqrt(2mE), where E are the obtained energy eigenvalues \n $l = {0},\, j = {1},\, N = {2},\, k\_max = {3},\, V_0 = {4}$ MeV, nr of contours$ = {5}$ \n The contour goes from Im[k]=-{6} to Im[k]=-{7}".format(l, j, order, k_max, problem.V0, cont_nrs, peaks_y[0], peaks_y[-1]), fontsize=20)
plt.xlabel('Re[k]', fontsize=20)
plt.ylabel('Im[k]', fontsize=20)
print 'V0 =', V0,', N =', order,', nr of contours =', cont_nrs
k_res=sp.empty(cont_nrs, 'complex')
#eigvecs_all=sp.empty((order+2, order+2, cont_nrs), 'complex')
contours = sp.empty((order, cont_nrs), 'complex')
k = sp.empty((order, cont_nrs), 'complex')
#eigvecs_res=sp.empty((order+2, cont_nrs), 'complex')
#absq_wavefunctions_res=sp.empty((2000, cont_nrs))
#eigvecs_res=sp.empty(order, cont_nrs, 'complex')
for m, peak_y in enumerate(peaks_y):
    print cont_nrs-m, 'to go'
    contour = calc.triangle_contour(peak_x, peak_y, k_max, order/3)
    ks, _ = contour
    H = calc.contour_hamiltonian(mom.H_element_contour, contour, args)
    eigvals, eigvecs = calc.energies(H)
    res = res_index(eigvecs)
    k_res[m] = sp.sqrt(2*problem.mass*eigvals[res])
    k[:,m]=sp.sqrt(2*problem.mass*eigvals)
    contours[:,m]=ks
    #absq_wavefunctions_res[:,m] = r**2 * absq(reswf(r))
    #eigvecs_all[:,:,m]=eigvecs
    #eigvecs_res[:,m] = eigvecs[:,res]

#k=sp.sqrt(2*problem.mass*eigvals)
plt.plot(sp.real(k), sp.imag(k), 'ko', markersize=3)
plt.plot(sp.real(ks), sp.imag(ks), 'ro', markersize=3)
plt.plot(sp.real(k_res), sp.imag(k_res), 'go')

#plt.axis([0.16, 0.18, -0.038, -0.028])
#file_name_2 = './mega__sweep_47_zoom_2.pdf'
#plt.savefig(file_name_2, format="pdf" )

#plt.axis([0, k_max, -peaks_y[-1], 0])
#file_name = "./mega_sweep_47_2.pdf"
#plt.savefig(file_name, format="pdf" )

print 'done'
#res=res_index(eigvecs)

#plt.figure(2)
#plt.xlabel('Re[k]', fontsize=20)
#plt.ylabel('$|\Phi_k|$', fontsize=20)
#plt.title('Eigenvectors in k-space \n N = {0}, V_0 = {1} MeV, peak Im[k]=-{2}'.format(order, V0, peaks_y[-1]), fontsize=20)
#plt.plot(sp.real(ks), abs(eigvecs))
#plt.plot(sp.real(ks), abs(eigvecs[:,res]), 'k', linewidth=3)

#rmax=100
#r = sp.linspace(1e-1, rmax, 2000)

#basis_function = mom.gen_basis_function(problem, l=l, j=j)

#reswf = calc.gen_wavefunction(eigvecs[:,res], basis_function, contour)
#reswf = calc.normalize(reswf, 0, 20, weight= lambda r: r**2)

#plt.figure(3)
#plt.title('Coordinate space wavefunction \n $N = {0},\, V_0 = {1}$ MeV, peak Im$[k] = {2}$'.format(order, V0, -peaks_y[-1]), fontsize=20)
#plt.xlabel(r'$r$ [fm]', fontsize=20)
#plt.ylabel(r'$r^2|\Psi(r)|^2$', fontsize=20)

#plt.plot(r, r**2 * absq(reswf(r)), 'k', linewidth=2)


plt.show()
