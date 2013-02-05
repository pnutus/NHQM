from math import *
import scipy as sp
import numpy as np
from scipy.fftpack import fft, ifft

from numpy import sin, linspace, pi
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, arange



@np.vectorize

def testFunc( x, y ):
    return (x % 2) == (y % 2)
    
n = 4
matris = sp.mat( sp.empty( (n, n) ))  
for x in xrange( n ):
    for y in xrange( n ):
        matris[ x, y ] = testFunc( x, y )
        
transM = sp.fftpack.fft2( matris )
h = []
for z in xrange(10):
   h.append( z )
   


def plotSpectrum(y,Fs):
 """
 Plots a Single-Sided Amplitude Spectrum of y(t)
 """
 n = len(y) # length of the signal
 k = arange(n)
 T = n/Fs
 frq = k/T # two sides frequency range
 frq = frq[range(n/2)] # one side frequency range

 Y = fft(y)/n # fft computing and normalization
 print Y
 Y = Y[range(n/2)]
 
 plot(frq,abs(Y),'r') # plotting the spectrum
 xlabel('Freq (Hz)')
 ylabel('|Y(freq)|')

Fs = 150.0;  # sampling rate
Ts = 1.0/Fs; # sampling interval
t = arange(0,1,Ts) # time vector

ff = 5;   # frequency of the signal
y = sin(2*pi*ff*t)
x = 0.7*sin(2*pi*3*t)

for i in xrange(len(y)):
    y[i] += x[i]

subplot(2,1,1)
plot(t,y)
xlabel('Time')
ylabel('Amplitude')
subplot(2,1,2)
plotSpectrum(y,Fs)
show()
            