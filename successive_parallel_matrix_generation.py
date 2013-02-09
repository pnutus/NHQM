import numpy as np
import scipy as sp
import time
from numpy import *
import parallel_matrix as pm
import os

def find_existing_matrix( matrix_file_name_basis, matrix_file_path, matrix_file_format,  sought_dim = -1 , verbose = False ):

    """Scans a folder and find a matrix stored in a file with a given name stem
    and returns the most suitable matrix, if no n is specified the largest is returned
    
    returns a matrix and its file path (NULL if no suitable matrix was found)"""
    
    output_strings = []
    
    #Redundancy check for directory names, making sure it ends with the required '/'
    if matrix_file_path[ len( matrix_file_path )-1 ] != '/':
        matrix_file_path = matrix_file_path + '/'
    
    existing_matrix_dimension = 0 
    temporary_dimension = -1
    
    #Scanning the directory and analyzing the sizes of the matrices, if possible chosing a matrix 
    #just larger than n, or if no such matrix exist returning the largest matrix
    os.chdir( matrix_file_path ) 
    for files in os.listdir( "." ):
        if files.rfind( matrix_file_name_basis ) == 0: #checks whether the filename starts with desired substring
            output_strings.append( 'applicable file: ' + files )
            temporary_dimension = int( files[ len( matrix_file_name_basis ):( len(files) - 4 )] )
            
            if temporary_dimension > existing_matrix_dimension:
                existing_matrix_dimension = temporary_dimension
                if existing_matrix_dimension >= sought_dim and sought_dim > 0:
                    break
        
    #loads the matrix from file    
    if existing_matrix_dimension > 0:            
        existing_matrix_path = matrix_file_path + matrix_file_name_basis + str( existing_matrix_dimension ) + matrix_file_format 
        output_strings.append( 'loading: ' + existing_matrix_path ) 
        result = sp.mat( sp.genfromtxt( existing_matrix_path ))
    else:
        result = sp.mat( sp.empty( (0,0) ))
        existing_matrix_path = 'NULL'
        
    if verbose:
        for arg in output_strings:
            print arg    
    
    return result, existing_matrix_path

"""def generate_succ_para_matrix_from_matrix( matrix_dim, element_generating_function, existing_matrix, matrix_file_name_basis, matrix_file_path, matrix_file_format, verbose = False):
    
    #Generates a new matrix using some function of the indices and copying existing values from a 
    #provided matrix. Saves the result in a txt file according the the matrixFile... parameters
    
    #Returns a matrix of desired dimensions
   
    #Redundancy check for directory names, making sure it ends with the required '/'
    if matrix_file_path[ len( matrix_file_path ) - 1 ] != '/':
        matrix_file_path = matrix_file_path + '/'
    
    existing_matrix_dimension,_ = existing_matrix.shape
    
    matris = parallel_matrix.parallel_matrix( element_generating_function, matrix_dim, False, existing_matrix_dimension )
        
    output_strings = []
    
    #Copies existing values from the existing matrix so as to avoid having to calculate the values
    #in a costly manner
    if existing_matrix_dimension > 0: 
        for i in xrange( min( existing_matrix_dimension, matrix_dim ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existing_matrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    #If the existing matrix is smaller than the required one the remaining valuues are calculated
    #using provided function.            
    if matrix_dim > existing_matrix_dimension:
    
    #saves the matrix to a .txt file
    save_matrix_path = matrix_file_path + matrix_file_name_basis + str(matrix_dim) + matrix_file_format
    output_strings.append( 'saving to: ' + save_matrix_path )
    savetxt(save_matrix_path, matris)
    
    if verbose:
        for arg in output_strings:
            print arg
    
    return matris"""
    
def generate_succ_para_matrix( matrix_dim, element_generating_function, matrix_file_name_basis, matrix_file_path, matrix_file_format, verbose = False):

    """Generates a new matrix using some function of the indices and if possible copies existing values from a matrix saved in a .txt file. Saves the result in a .txt file according the the matrixFile... parameters
    
    Returns a matrix of desired dimensions"""
   
    #Redundancy check for directory names, making sure it ends with the required '/'   
    if matrix_file_path[ len( matrix_file_path ) - 1 ] != '/':
        matrix_file_path = matrix_file_path + '/'
        
    t = time.time()
    existing_matrix, old_matrix_path = find_existing_matrix( matrix_file_name_basis, matrix_file_path, matrix_file_format,  matrix_dim, verbose ) 
    print "loading existing matrix:"   
    print time.time() - t 
    
    existing_matrix_dimension,_ = existing_matrix.shape
    
    t = time.time()
    matris = pm.parallel_matrix( element_generating_function, matrix_dim, False, existing_matrix_dimension)
    print "creating new matrix elements in parallel:"
    print time.time() - t
    
    output_strings = []
    
    #Copies existing values from the existing matrix so as to avoid having to calculate the values
    #in a costly manner
    t = time.time()
    if existing_matrix_dimension > 0: 
        for i in xrange( min( existing_matrix_dimension, matrix_dim ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existing_matrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
                        
    print "copying existing elements:"
    print time.time() - t                    
                    
    save_matrix_path = matrix_file_path + matrix_file_name_basis + str(matrix_dim) + matrix_file_format
    #verify that the new and old matrices are not the same to avoid needless file operations
    #and make sure to only keep the file with the most information
    if ( save_matrix_path != old_matrix_path and matrix_dim > existing_matrix_dimension ):                
        #saves the matrix to a .txt file
        output_strings.append( 'saving to: ' + save_matrix_path )
        savetxt(save_matrix_path, matris) 
    
        """might want to double check whether the save was succesful before deleting 
        the old matrix file"""
    
        #deletes the old matrix file
        if old_matrix_path != 'NULL':
            output_strings.append( 'deleting: ' + old_matrix_path )    
            os.remove(old_matrix_path)
    
    
    if verbose:
        for arg in output_strings:
            print arg
    
    return matris    
