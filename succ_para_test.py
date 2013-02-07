import time
import find_existing_matrix
import He5
import HO_basis
import successive_parallel_matrix_generation
import successive_matrix_generation
import parallel_matrix
import scipy as sp



matris_path = '/Users/Spill/Documents/kandidatarbete/matrisdata'
matris_format = '.txt'
matris_namn = 'He5_matris_para_test'

def one(i,j):
    return 1

#element generating function
def H_element_wrapper(n,nprim):
    #lement(n, n_prim, l, s, omega, V):
    return HO_basis.H_element(n,nprim,0,0.5,1,He5.V) + n
    
def seq_mat(f, dim):
    mat = sp.mat ( sp. empty((dim,dim)) )
    for r in xrange(dim):
        for k in xrange(dim):
            mat[r,k] = f(r,k)
    
    return mat    

dim = 50
"""
x = time.time()

mat = successive_parallel_matrix_generation.generate_succ_para_matrix( dim, H_element_wrapper, matris_namn, matris_path, matris_format )

print "parallel: " + str( time.time() - x )


matris_namn = 'He5_matris_test'

y = time.time()

hamilton = successive_matrix_generation.generate_succ_matrix( dim, H_element_wrapper, matris_namn, matris_path, matris_format, False )

print "sekvens: " + str( time.time() - y )

print mat - hamilton


"""
t = time.time()    
#sekv = seq_mat(H_element_wrapper, dim)
sekv_t = t - time.time() 
print sekv_t

t = time.time()
para = parallel_matrix.parallel_matrix(H_element_wrapper, dim)
para_t = t - time.time()    
print para_t

print sekv - para

t = time.time()
parax = parallel_matrix.parallel_matrix(H_element_wrapper, dim, False, dim // 2)
parax_t = t - time.time()    
print parax_t

#print sekv - para

print "sekventiellt:"
print sekv_t
print "parallellt:"
print para_t   
print "paraex:"
print parax_t
            
            

