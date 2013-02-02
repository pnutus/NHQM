from math import *

import numpy as np
import scipy as sp
from numpy import *

import os

def findExistingMatrix(n,matrixFileNameBasis,matrixFilePath,matrixFileFormat)
    
    existingMatrixDimension = 0 
    
    os.chdir(matrixFilePath) 
    for files in os.listdir("."):
        temporaryDimension = int( files[len(matrixFileNameBasis):(len(files)-4)] )
            
    if temporaryDimension <= n:
        existingMatrixDimension = temporaryDimension
    
    
    existingMatrixPath = matrixFilePath + matrixFileNameBasis + str(existingMatrixDimension) + matrixFileFormat 
    
    if existingMatrixDimension != 0:      
        return genfromtxt(existingMatrixPath)
    else:
        return 0
