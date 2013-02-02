from math import *

import numpy as np
import scipy as sp
from numpy import *

import os

def generateSuccMatrix(n,generatingFunction,existingMatrix,matrixFileNameBasis,matrixFilePath,matrixFileFormat)    
    matris = sp.zeros((n,n))
    
    if existingMatrixDimension > 1:    
        for i in range(0,existingMatrixDimension):
                for j in range(0,i+1):
                    matris[i,j] = existingMatrix[i,j]
                    if not (i == j):
                        matris[j,i]=matris[i,j]
    
    if existingMatrixDimension == 1:
        matris[1,1] = existingMatrix
    
    for i in range(existingMatrixDimension, n):
        for j in range(0 , existingMatrixDimension):
            matris[j,i] = matris[i,j] = elementGeneratingFunction(i,j)

    
    for i in range(existingMatrixDimension , n):        
        for j in range(existingMatrixDimension, i + 1): 
            matris[i,j] = elementGeneratingFunction(i,j)       
            if ( i != j):
                matris[j,i] = matris[i,j]
    
    saveMatrixPath = matrixFilePath + matrixFileNameBasis + str(n) + matrixFileFormat
    savetxt(saveMatrixPath, matris)
    
    
    return matris

