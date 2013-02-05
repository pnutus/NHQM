from math import *

import numpy as np
import scipy as sp
from numpy import *

import os

def findExistingMatrix( matrixFileNameBasis, matrixFilePath, matrixFileFormat,  n = -1 , verbose = False ):

    """Scans a folder and find a matrix stored in a file with a given name stem
    and returns the most suitable matrix, if no n is specified the largest is returned
    
    returns a matrix and its file path (NULL if no suitable matrix was found)"""
    
    outputStrings = []
    
    #Redundancy check for directory names, making sure it ends with the required '/'
    if matrixFilePath[ len( matrixFilePath )-1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
    
    existingMatrixDimension = 0 
    temporaryDimension = -1
    
    #Scanning the directory and analyzing the sizes of the matrices, if possible chosing a matrix 
    #just larger than n, or if no such matrix exist returning the largest matrix
    os.chdir( matrixFilePath ) 
    for files in os.listdir( "." ):
        if files.rfind( matrixFileNameBasis ) == 0: #checks whether the filename starts with desired substring
            outputStrings.append( 'applicable file: ' + files )
            temporaryDimension = int( files[ len( matrixFileNameBasis ):( len(files) - 4 )] )
            
            if temporaryDimension > existingMatrixDimension:
                existingMatrixDimension = temporaryDimension
                if existingMatrixDimension >= n and n > 0:
                    break
        
    #loads the matrix from file    
    if existingMatrixDimension > 0:            
        existingMatrixPath = matrixFilePath + matrixFileNameBasis + str( existingMatrixDimension ) + matrixFileFormat 
        outputStrings.append( 'loading: ' + existingMatrixPath ) 
        result = sp.mat( sp.genfromtxt( existingMatrixPath ))
    else:
        result = sp.mat( sp.empty( (0,0) ))
        existingMatrixPath = 'NULL'
        
    if verbose:
        for arg in outputStrings:
            print arg    
    
    return result, existingMatrixPath

def generateSuccMatrixFromMatrix( n, elementGeneratingFunction, existingMatrix, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose = False):
    
    """Generates a new matrix using some function of the indices and copying existing values from a 
    provided matrix. Saves the result in a txt file according the the matrixFile... parameters
    
    Returns a matrix of desired dimensions"""
   
    #Redundancy check for directory names, making sure it ends with the required '/'
    if matrixFilePath[ len( matrixFilePath ) - 1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
    
    existingMatrixDimension,_ = existingMatrix.shape
    
    matris = sp.mat(sp.empty( (n,n) ))
        
    outputStrings = []
    
    #Copies existing values from the existing matrix so as to avoid having to calculate the values
    #in a costly manner
    if existingMatrixDimension > 0: 
        for i in xrange( min( existingMatrixDimension, n ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existingMatrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    #If the existing matrix is smaller than the required one the remaining valuues are calculated
    #using provided function.            
    if n > existingMatrixDimension:
            
        for i in xrange( existingMatrixDimension, n ):
            for j in xrange( 0 , existingMatrixDimension ):
                outputStrings.append( 'calculating element: ' + str(i) + ', ' + str(j) )
                matris[j,i] = matris[i,j] = elementGeneratingFunction( i, j )
        
        
        for i in xrange( existingMatrixDimension , n ):        
            for j in xrange( existingMatrixDimension, i + 1 ): 
                outputStrings.append( 'calculating element: ' + str(i) + ', ' + str(j) )
                matris[i,j] = elementGeneratingFunction(i,j)       
                if ( i != j):
                    matris[j,i] = matris[i,j]
    
    #saves the matrix to a .txt file
    saveMatrixPath = matrixFilePath + matrixFileNameBasis + str(n) + matrixFileFormat
    outputStrings.append( 'saving to: ' + saveMatrixPath )
    savetxt(saveMatrixPath, matris)
    
    if verbose:
        for arg in outputStrings:
            print arg
    
    return matris
    
def generateSuccMatrix( n, elementGeneratingFunction, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose = False):

    """Generates a new matrix using some function of the indices and if possible copies existing values from a matrix saved in a .txt file. Saves the result in a .txt file according the the matrixFile... parameters
    
    Returns a matrix of desired dimensions"""
   
    #Redundancy check for directory names, making sure it ends with the required '/'   
    if matrixFilePath[ len( matrixFilePath ) - 1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
        
    existingMatrix, oldMatrixPath = findExistingMatrix( matrixFileNameBasis, matrixFilePath, matrixFileFormat,  n, verbose )    
    
    existingMatrixDimension,_ = existingMatrix.shape
    
    matris = sp.mat(sp.empty( (n,n) ))
        
    outputStrings = []
    
    #Copies existing values from the existing matrix so as to avoid having to calculate the values
    #in a costly manner
    if existingMatrixDimension > 0: 
        for i in xrange( min( existingMatrixDimension, n ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existingMatrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    #If the existing matrix is smaller than the required one the remaining valuues are calculated
    #using provided function.            
    if n > existingMatrixDimension:
            
        for i in xrange( existingMatrixDimension, n ):
            for j in xrange( 0 , existingMatrixDimension ):
                outputStrings.append( 'calculating element: ' + str(i) + ', ' + str(j) )
                matris[j,i] = matris[i,j] = elementGeneratingFunction( i, j )
        
        
        for i in xrange( existingMatrixDimension , n ):        
            for j in xrange( existingMatrixDimension, i + 1 ): 
                outputStrings.append( 'calculating element: ' + str(i) + ', ' + str(j) )
                matris[i,j] = elementGeneratingFunction(i,j)       
                if ( i != j):
                    matris[j,i] = matris[i,j]
                    
                    
    saveMatrixPath = matrixFilePath + matrixFileNameBasis + str(n) + matrixFileFormat
    #verify that the new and old matrices are not the same to avoid needless file operations
    #and make sure to only keep the file with the most information
    if ( saveMatrixPath != oldMatrixPath and n > existingMatrixDimension ):                
        #saves the matrix to a .txt file
        outputStrings.append( 'saving to: ' + saveMatrixPath )
        savetxt(saveMatrixPath, matris) 
    
        """might want to double check whether the save was succesful before deleting 
        the old matrix file"""
    
        #deletes the old matrix file
        if oldMatrixPath != 'NULL':
            outputStrings.append( 'deleting: ' + oldMatrixPath )    
            os.remove(oldMatrixPath)
    
    
    if verbose:
        for arg in outputStrings:
            print arg
    
    return matris    