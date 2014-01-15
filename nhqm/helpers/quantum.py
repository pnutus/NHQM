# This Python file uses the following encoding: utf-8
from __future__ import division
import scipy as sp
from scipy import linalg, integrate, floor, sqrt, arange
from scipy.special.specfun import csphjy
from scipy.misc import factorial
from nhqm.helpers.decorators import memoize


def j_l(l, x):
    """
    The Spherical Bessel function j_l.
    
    Reference:
    http://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn.2C_yn
    """
    L = 1 if l == 0 else l
    _, j_l, _, _, _ = csphjy(L, x)
    return j_l[l]
    
def absq(x):
    """
    Computes the absolute value squared. Faster than norm(x)^2 since
    it avoids the square root operation.
    """
    return sp.real(x)**2 + sp.imag(x)**2

def normalize(x, norm=sp.norm):
    """
    Normalizes an array using supplied norm.
    """
    return x / sp.norm(x)

def L2_norm(f, a, b, weight = lambda x: 1):
    """
    Computes the L2 norm of a function f given limits a and b 
    and a weight for the L2 space.
    """
    def integrand(x):
        return absq(f(x)) * weight(x);
    (N, _) = sp.integrate.quad(integrand, a, b)
    return N

def L2_normalize(f, a, b, weight = lambda x: 1):
    """
    Normalizes a function over an L2 space. See L2_norm. 
    """
    N = L2_norm(f, a, b, weight = weight)
    return lambda x: f(x) / sp.sqrt(N)

#
# Angular momentum stuff
#

# Limit on number of iterations for computing C-G coefficients
CG_LIMIT = 50

@memoize
def clebsch_gordan(j1, j2, m1, m2, J, M):
    """
    Computes the Clebsch-Gordan coefficient
    <j1 j2; m1 m2|j1 j2; J M>.

    For reference see
    http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients.
    """
    if M != m1 + m2 or not abs(j1 - j2) <= J <= j1 + j2:
        return 0
    c1 = sp.sqrt((2*J+1) * factorial(J+j1-j2) * factorial(J-j1+j2) * \
              factorial(j1+j2-J)/factorial(j1+j2+J+1))
    c2 = sp.sqrt(factorial(J+M) * factorial(J-M) * factorial(j1-m1) * \
              factorial(j1+m1) * factorial(j2-m2) * factorial(j2+m2))
    c3 = 0.
    for k in range(CG_LIMIT):
        use = True
        d = [0, 0, 0, 0, 0]
        d[0] = j1 + j2 - J - k
        d[1] = j1 - m1 - k
        d[2] = j2 + m2 - k
        d[3] = J - j2 + m1 + k
        d[4] = J - j1 -m2 + k
        prod = factorial(k)
        for arg in d:
            if arg < 0:
                use = False
                break
            prod *= factorial(arg)
        if use:
            c3 += (-1)**k/prod
    return c1*c2*c3
    

def wigner3j(j1,j2,j3,m1,m2,m3):
#======================================================================
# Wigner3j.m by David Terr, Raytheon, 6-17-04
#
# Compute the Wigner 3j symbol using the Racah formula [1]. 
#
# Usage: 
# from wigner import Wigner3j
# wigner = Wigner3j(j1,j2,j3,m1,m2,m3)
#
#  / j1 j2 j3 \
#  |          |  
#  \ m1 m2 m3 /
#
# Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: 
# http://mathworld.wolfram.com/Wigner3j-Symbol.html
#======================================================================

    # Error checking
    if ( ( 2*j1 != floor(2*j1) ) | ( 2*j2 != floor(2*j2) ) | ( 2*j3 != floor(2*j3) ) | ( 2*m1 != floor(2*m1) ) | ( 2*m2 != floor(2*m2) ) | ( 2*m3 != floor(2*m3) ) ):
        raise ValueError('All arguments must be integers or half-integers.')
        

    # Additional check if the sum of the second row equals zero
    if (   m1 + m2 + m3 != 0
        or j1 - m1 != floor( j1 - m1 )
        or j2 - m2 != floor( j2 - m2 )
        or j3 - m3 != floor( j3 - m3 )
        or j3 - m3 != floor( j3 - m3 )
        or j3 > j1 + j2 or j3 < abs(j1 - j2)
        or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3
        ):
        return 0

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max( 0, max( t1, t2 ) )
    tmax = min( t3, min( t4, t5 ) )
    tvec = arange(tmin, tmax+1, 1)

    wigner = 0

    for t in tvec:
        wigner += (-1)**t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )

    return wigner * (-1)**(j1-j2-m3) * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1) * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) )


def wigner6j(j1,j2,j3,J1,J2,J3):
#======================================================================
# Calculating the Wigner6j-Symbols using the Racah-Formula                
# Author: Ulrich Krohn                                            
# Date: 13th November 2009
#                                                                         
# Based upon Wigner3j.m from David Terr, Raytheon                         
# Reference: http://mathworld.wolfram.com/Wigner6j-Symbol.html            
#
# Usage: 
# from wigner import Wigner6j
# WignerReturn = Wigner6j(j1,j2,j3,J1,J2,J3)
#
#  / j1 j2 j3 \
# <            >  
#  \ J1 J2 J3 /
#
#======================================================================

    # Check that the js and Js are only integer or half integer
    if ( ( 2*j1 != round(2*j1) ) | ( 2*j2 != round(2*j2) ) | ( 2*j2 != round(2*j2) ) | ( 2*J1 != round(2*J1) ) | ( 2*J2 != round(2*J2) ) | ( 2*J3 != round(2*J3) ) ):
        print 'All arguments must be integers or half-integers.'
        return -1
    
# Check if the 4 triads ( (j1 j2 j3), (j1 J2 J3), (J1 j2 J3), (J1 J2 j3) ) satisfy the triangular inequalities
    if ( ( abs(j1-j2) > j3 ) | ( j1+j2 < j3 ) | ( abs(j1-J2) > J3 ) | ( j1+J2 < J3 ) | ( abs(J1-j2) > J3 ) | ( J1+j2 < J3 ) | ( abs(J1-J2) > j3 ) | ( J1+J2 < j3 ) ):
        print '6j-Symbol is not triangular!'
        return 0
    
    # Check if the sum of the elements of each traid is an integer
    if ( ( 2*(j1+j2+j3) != round(2*(j1+j2+j3)) ) | ( 2*(j1+J2+J3) != round(2*(j1+J2+J3)) ) | ( 2*(J1+j2+J3) != round(2*(J1+j2+J3)) ) | ( 2*(J1+J2+j3) != round(2*(J1+J2+j3)) ) ):
        print '6j-Symbol is not triangular!'
        return 0
    
    # Arguments for the factorials
    t1 = j1+j2+j3
    t2 = j1+J2+J3
    t3 = J1+j2+J3
    t4 = J1+J2+j3
    t5 = j1+j2+J1+J2
    t6 = j2+j3+J2+J3
    t7 = j1+j3+J1+J3

    # Finding summation borders
    tmin = max(0, max(t1, max(t2, max(t3,t4))))
    tmax = min(t5, min(t6,t7))
    tvec = arange(tmin,tmax+1,1)
        
    # Calculation the sum part of the 6j-Symbol
    WignerReturn = 0
    for t in tvec:
        WignerReturn += (-1)**t*factorial(t+1)/( factorial(t-t1)*factorial(t-t2)*factorial(t-t3)*factorial(t-t4)*factorial(t5-t)*factorial(t6-t)*factorial(t7-t) )

    # Calculation of the 6j-Symbol
    return WignerReturn*sqrt( TriaCoeff(j1,j2,j3)*TriaCoeff(j1,J2,J3)*TriaCoeff(J1,j2,J3)*TriaCoeff(J1,J2,j3) )


def TriaCoeff(a,b,c):
    # Calculating the triangle coefficient
    return factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1))
