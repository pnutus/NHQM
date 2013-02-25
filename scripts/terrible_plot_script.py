 
import time
import numpy as np

problem = He5.problem
basis = mom
#print np.sqrt(2-1j)

def find(v,x):
    for i in xrange(len(v)):
        if (v[i] == x):
           # print i
            return i
    #print str(-1)        
    return -1


                 
def kontroll_samma_kontur(order, k_max, peak_pos, peak_amp, zero_point, l=0, j=0.5):
    step_size = k_max / order #only step_size = 1 works with the current contour implementation...

    r = sp.linspace(0, 10, 200)

    kontur_a = berg.gen_berggren_contour( peak_pos, peak_amp, k_max, order, zero_point )
    kontur_b = berg.gen_berggren_contour( peak_pos, peak_amp, k_max, order, zero_point )

    a = 0
    kont_hamilton = sp.empty( (order,order), complex ) # complex elements?
    for p in kontur_a:
        b = 0
        for p_prim in kontur_b:
            kont_hamilton[a,b] = berg.H_complex_element(p, p_prim, He5.problem, step_size, l, j)
            b += 1
        a += 1
        kontur_b = berg.gen_berggren_contour( peak_pos, peak_amp, k_max, order, zero_point )  
    
    
    #H = calc.hamiltonian(basis.H_element, args=(problem, step_size, l, j), order=order)
    H=[[0]]

    return kont_hamilton , H

def plott_kontur( peak_position, peak_neg_amplitude, k_max, points, first_zero_position = -1  ):
    k = berg.gen_berggren_contour( peak_position, peak_neg_amplitude, k_max, points, first_zero_position )
    res = []
    for i in k:
         res.append(i)
         
    re = sp.real(res)
    im = sp.imag(res)
    
    plt.xlabel('Re')
    plt.ylabel('Im')
    #plt.title('Berggren kontur')
    

    plt.plot(re, im)

    lin = sp.linspace(0,k_max)
    noll = sp.zeros((len(lin),1))
    
    plt.plot(lin, noll, ls='-')
    plt.show()
    
    
if __name__ == '__main__':
    
    problem.V0=-47.0
    
    k_max = 2.5
    pp = .185
    #pa = 0.0001
    zp = 0.4
    l = 1
    j = 1.5
    ###order is the size of the matrix (number of mesh points)
    order = 100
    ###porder is the number of contours that is included in the ''sweep''
    porder=1
    ###palist contains all the ''peak amplitudes'' of the contour (that is, how far down the contour goes)
    ###note that there is currently a bug such that you cannot use pa=0. Instead you have to start from, say, pa=0.000001  
    palist=sp.linspace(0.05, 0.05, porder)
    
    Re_energies=sp.empty((porder,order))
    Im_energies=sp.empty((porder,order))
    Re_k=sp.empty((porder,order))
    Im_k=sp.empty((porder,order))
    for k, pa in enumerate(palist):
        print porder-k, 'to go'
        h_c, h =kontroll_samma_kontur(order, k_max, pp, pa, zp, l, j)
        eigval_c, eigvec_c = calc.energies(h_c)
        k_c = np.sqrt(eigval_c * 2 * problem.mass)
        Re_k[k,:]=sp.real(k_c)
        Im_k[k,:]=sp.imag(k_c)
        Re_energies[k,:]=sp.real(eigval_c)
        Im_energies[k,:]=sp.imag(eigval_c)
       
    print 'done'

    
    plt.figure(1)
    plt.xlabel('Re[E]', fontsize=20)
    plt.ylabel('Im[E]', fontsize=20)
    for m in range(0,order):
        plt.plot(Re_energies[:,m],Im_energies[:,m],'ko')
        
    plt.figure(2)
    plt.title("The figure shows the motion of the poles k=sqrt(2mE),\n where E are the obtained energy eigenvalues,\n as the contour is continuously varied \n $l = {0},\, j = {1},\, N = {2},\, kmax = {3},\, V_0 = {4}$ MeV, nr of contours$ = {5}$".format(l, j, order, k_max, problem.V0, porder), fontsize=20)
    plt.xlabel('Re[k]', fontsize=20)
    plt.ylabel('Im[k]', fontsize=20)
    for m in range(0,order):
        plt.plot(Re_k[:,m],Im_k[:,m],'ko', markersize=4)

    plt.show()
   
