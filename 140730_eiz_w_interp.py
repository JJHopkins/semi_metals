#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show
import pyreport
#from scipyx.interpolate import interp1d

x_eV,y_eV = np.loadtxt('data/water-LDS-02.csv',delimiter = ',',unpack=True, usecols = [0,1])
x_eV,y_eV = np.loadtxt('data/AL2O3-VUV-L.csv',delimiter = ',',unpack=True, usecols = [0,1])
x = x_eV/0.159
y = y_eV
#x = x_eV[1:]/0.159
#y = y_eV[1:]
f = np.zeros(shape = (len(x)))
#f = np.zeros(shape = (len(x),len(y)))

X = np.arange(0,12501)

for i,ex in enumerate(X):
    f = np.interp(X,x,y)
    print i#y(X)
    #f[i] = interp1d(x, y, kind = 'cubic')

pl.figure()
pl.loglog(x,y,'b:+',X,f,'g--x')
pl.show()

np.savetxt( "data/140730_eiz_w.txt", f[:500])
