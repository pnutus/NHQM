import time
import He5
import HO_basis as HO
import successive_parallel_matrix_generation as spmg
import successive_matrix_generation as smg
import parallel_matrix as pm
import scipy as sp



matris_path = '/Users/Spill/Documents/kandidatarbete/matrisdata'
matris_format = '.txt'
matris_namn = 'He5_matris_para_test'

def one(i,j):
    return 1

#element generating function
def H_element_wrapper(n,nprim):
    #lement(n, n_prim, l, s, omega, V):
    return HO.H_element(n,nprim,0,0.5,1,He5.V) + n
        
def seq_mat(f, dim):
    mat = sp.mat ( sp. empty((dim,dim)) )
    for r in xrange(dim):
        for k in xrange(dim):
            mat[r,k] = f(r,k)
    
    return mat    

dim = 100
"""
x = time.time()

mat = spmg.generate_succ_para_matrix( dim, H_element_wrapper, matris_namn, matris_path, matris_format )

print "parallel: " + str( time.time() - x )


matris_namn = 'He5_matris_test'

y = time.time()

hamilton = smg.generate_succ_matrix( dim, H_element_wrapper, matris_namn, matris_path, matris_format, False )

print "sekvens: " + str( time.time() - y )

print mat - hamilton



t = time.time()
para = pm.parallel_matrix(H_element_wrapper, dim)
para_t = time.time() - t   
print para_t

t = time.time()    
sekv = seq_mat(H_element_wrapper, dim)
sekv_t = time.time() - t
print sekv_t

print sekv - para

t = time.time()
parax = pm.parallel_matrix(H_element_wrapper, dim, False, dim -1)
parax_t = time.time() - t   
print parax_t

#print sekv - para

print "sekventiellt:"
print sekv_t
print "parallellt:"
print para_t   
print "paraex:"
print parax_t
            
t = time.time()    
sekv2 = seq_mat(H_element_wrapper, dim)
sekv_t2 = time.time() - t

t = time.time()
para2 = pm.parallel_matrix(H_element_wrapper, dim)
para_t2 = time.time() - t   

t = time.time()
parax2 = pm.parallel_matrix(H_element_wrapper, dim, False, dim -1)
parax_t2 = time.time() - t   


print "sekventiellt: 2"
print sekv_t2
print "parallellt: 2"
print para_t2  
print "paraex: 2"
print parax_t2            """
            
            

t = time.time()
mat = pm.parallel_matrix(H_element_wrapper, dim)     
print time.time() - t      
     
t = time.time()
#mat = seq_mat(H_element_wrapper, dim)     
print time.time() - t      
     
     

