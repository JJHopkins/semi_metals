#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import *
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

Ls = np.arange(1e-9,450e-9,1e-9)  # separation distance between 2 cyclinders
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]


A65  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_65_sum.txt') # LDS in perpendicular direction
A91  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_91_sum.txt') # LDS in perpendicular direction
A93  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_93_sum.txt') # LDS in perpendicular direction
A290 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_290_sum.txt') # LDS in perpendicular direction 

A265  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_65_sum.txt') # LDS in perpendicular direction
A291  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_91_sum.txt') # LDS in perpendicular direction
A293  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_93_sum.txt') # LDS in perpendicular direction
A2290=  np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_290_sum.txt') # LDS in perpendicular direction 

An65   =  np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_65_n.txt') # LDS in perpendicular direction
An91   =  np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_91_n.txt') # LDS in perpendicular direction
An93   =  np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_93_n.txt') # LDS in perpendicular direction
An290  =  np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_290_n.txt') # LDS in perpendicular direction 

A2n65  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_65_n.txt') # LDS in perpendicular direction
A2n91  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_91_n.txt') # LDS in perpendicular direction
A2n93  = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_93_n.txt') # LDS in perpendicular direction
A2n290 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_290_n.txt') # LDS in perpendicular direction 

pl.figure()
#pl.plot(1e9*Ls,1e21*(A65 -(kbT/32)* An65[:,0]), 'k-',  label = r'[6,5] $\cal{A}^{(0)}$')
#pl.plot(1e9*Ls,1e21*(A265-(kbT/32)*A2n65[:,0]), 'k--', label = r'[6,5] $\cal{A}^{(2)}$')
#                                                                         
#pl.plot(1e9*Ls,1e21*( A91-(kbT/32)*An91  [:,0]),'r-',  label = r'[9,1] $\cal{A}^{(0)}$')
#pl.plot(1e9*Ls,1e21*(A291-(kbT/32)*A2n91 [:,0]),'r--', label = r'[9,1] $\cal{A}^{(2)}$')
#                                                                          
#pl.plot(1e9*Ls,1e21*( A290-(kbT/32)*An290 [:,0]),'g-', label = r'[29,0] $\cal{A}^{(0)}$')
#pl.plot(1e9*Ls,1e21*(A2290-(kbT/32)*A2n290[:,0]),'g--',label = r'[29,0] $\cal{A}^{(2)}$')

pl.plot(1e9*Ls,1e21* A290,'g-', label = r'[29,0] $\cal{A}^{(0)}$')
pl.plot(1e9*Ls,1e21*A2290,'g--',label = r'[29,0] $\cal{A}^{(2)}$')

pl.xlabel(r'$\ell\,\,\,[nm]$')
pl.ylabel(r'$\cal{A}^{(0,2)}\,\,\, [zJ]$')
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.show()

