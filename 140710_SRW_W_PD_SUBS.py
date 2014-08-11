#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb,trapz
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction

eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium

#eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction

#eiz_w = np.loadtxt('data/eiz_w_DD.txt') # LDS of water, intervening medium

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
#coeff = 0.159           # [eV]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]

# Matsubara frequencies
#ns = z/(0.159)#
ns = np.arange(0.,500.) 
zs = ns * coeff         

#Ls = np.arange(1e-9,450e-9,1e-9)  # separation distance between 2 cyclinders
#Ls = np.arange(1e-9,20e-9,0.1e-9)  # separation distance between 2 cyclinders
Ls = np.arange(1e-9,10e-9,1e-9)  # separation distance between 2 cyclinders

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta_par(par,med):
	return (par - med)/med

def Delta_prp(perp,med):
	return (perp - med)/(perp + med)

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

#p = np.zeros(shape = (len(Ls),len(ns)))
#A0 = np.zeros(shape = (len(Ls),len(ns)))
#A2 = np.zeros(shape = (len(Ls),len(ns)))
#
#a =  Aiz(eiz_x,eiz_z,eiz_w)
#delta_par = Delta_par(eiz_z,eiz_w)
#delta_prp = Delta_prp(eiz_x,eiz_w)
P = np.linspace(0.,10.)#Pn(eiz_w[j],zs[j],L)
f0 = np.zeros(shape = (len(p),len(T))) 
Ft0 = np.zeros(len(p)) 
# Calculate 1 \geq n \leq 500 terms
for i,P in enumerate(p):
    for j,t in enumerate(T):
        f0[i,j] =T[j]*np.exp(-2.*np.sqrt(p[i]*p[i]*T[j]*T[j]+p[i]*p[i]))/(T[j]*T[j]+1.)*(T[j]*T[j]+2.)(T[j]*T[j]+2.)
        Ft0 = np.trapz(f0, axis = 1)
    print Ft0

pl.figure()
pl.loglog(p,Ft0)
pl.show()

p = np.linspace(0.,10.,100)#Pn(eiz_w[j],zs[j],L)
t  = np.linspace(0.,20.,200)#, 1.+2.**17)
P,T = np.meshgrid(p,t, indexing = 'ij')
f = (T/(T**2+P**2))*np.exp(-2.*np.sqrt(T**2+P**2))*(T**2+2.*P**2)**2
F = romb(f)
print F

pl.figure()
pl.plot(P,F)
pl.show()

