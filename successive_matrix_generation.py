from math import *

import numpy as np
import scipy as sp
from numpy import *

import os

def findExistingMatrix( n, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose = 0 ):
    
    outputStrings = []
    if matrixFilePath[ len( matrixFilePath )-1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
    
    existingMatrixDimension = 0 
    temporaryDimension = -1
    
    os.chdir( matrixFilePath ) 
    for files in os.listdir( "." ):
        if files.rfind( matrixFileNameBasis ) == 0:
            outputStrings.append( 'applicable file: ' + files )
            temporaryDimension = int( files[ len( matrixFileNameBasis ):( len(files) - 4 )] )
            
            if temporaryDimension > existingMatrixDimension:
                existingMatrixDimension = temporaryDimension
                if existingMatrixDimension >= n:
                    break
            
    
    
        
    if existingMatrixDimension > 0:            
        existingMatrixPath = matrixFilePath + matrixFileNameBasis + str( existingMatrixDimension ) + matrixFileFormat 
        outputStrings.append( 'loading: ' + existingMatrixPath ) 
        result = sp.mat( sp.genfromtxt( existingMatrixPath ))
    else:
        result = sp.mat( sp.empty( (0,0) ))
        existingMatrixPath = 'NULL'
        
    if verbose:
        for arg in outputStrings:
            arg    
    
    return result, existingMatrixPath

def generateSuccMatrixFromMatrix( n, elementGeneratingFunction, existingMatrix, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose = 0):
   
    if matrixFilePath[ len( matrixFilePath ) - 1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
    
    existingMatrixDimension,_ = existingMatrix.shape
    
    matris = sp.mat(sp.empty( (n,n) ))
        
    outputStrings = []
    
    if existingMatrixDimension > 0: 
        for i in xrange( min( existingMatrixDimension, n ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existingMatrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    if existingMatrixDimension == 1:
            if existingMatrix == 0:
                matris[0,0] = elementGeneratingFunction( 1, 1 )
            else:
                matris[0,0] = existingMatrix
                
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
    outputStrings.append( 'saving to: ' + saveMatrixPath )
    savetxt(saveMatrixPath, matris)
    
    if verbose:
        for arg in outputStrings:
            print arg
    
    return matris
    
def generateSuccMatrix( n, elementGeneratingFunction, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose = 0):
   
   
    if matrixFilePath[ len( matrixFilePath ) - 1 ] != '/':
        matrixFilePath = matrixFilePath + '/'
        
    existingMatrix, oldMatrixPath = findExistingMatrix( n, matrixFileNameBasis, matrixFilePath, matrixFileFormat, verbose )    
    
    existingMatrixDimension,_ = existingMatrix.shape
    
    matris = sp.mat(sp.empty( (n,n) ))
        
    outputStrings = []
    
    if existingMatrixDimension > 0: 
        for i in xrange( min( existingMatrixDimension, n ) ):
                for j in xrange( 0, i + 1 ):
                    matris[i,j] = existingMatrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    if existingMatrixDimension == 1:
            if existingMatrix == 0:
                matris[0,0] = elementGeneratingFunction( 1, 1 )
            else:
                matris[0,0] = existingMatrix
                
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
    outputStrings.append( 'saving to: ' + saveMatrixPath )
    savetxt(saveMatrixPath, matris) 
    
    """might want to double check wheter the save was succesful before deleting the old matrix file"""
    if oldMatrixPath != 'NULL':
        outputStrings.append( 'deleting: ' + oldMatrixPath )    
        os.remove(oldMatrixPath)
    
    
    if verbose:
        for arg in outputStrings:
            print arg
    
    return matris    