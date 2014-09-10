#!/usr/bin/env python
import numpy as np
from scipy.integrate import romb

# Load LDS for 2 materials,medium
eiz_x_024 = np.loadtxt('data/eiz_x_024.txt') 
eiz_z_024 = np.loadtxt('data/eiz_z_024.txt') 
eiz_x_093 = np.loadtxt('data/eiz_x_093.txt') 
eiz_z_093 = np.loadtxt('data/eiz_z_093.txt') 
eiz_w = np.loadtxt('data/eiz_w.txt') 

# Values
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.)  # Matsubara sum
zs = ns * coeff         
r_1 = 1.6272e-9 # radius 24
r_2 = 0.4234e-9 # radius 93

thetas = np.arange(0.783,0.786,0.001) # mutual angle
Ls     = np.arange(1e-9,100e-9,1e-9)  # separation 

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Delta_perp(perp,med):
	return (perp - med)/(perp+med)

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(len(Ls))

# n = 0 terms
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a_1[0])*(1.+3.*a_2[0])
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a_1[0])*(1.-a_2[0])
Ft0_term0 = romb(f0_term0)
Ft2_term0 = romb(f2_term0)

for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a_1[j])*(1.+3.*a_2[j])*T*T*T*T\
                + 4.*(1.+2.*a_1[j]+2.*a_2[j]+3.*a_1[j]*a_2[j])*T*T \
                + 4.*(1.+a_1[j])*(1.+a_2[j]))
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a_1[j])*(1.-a_2[j]))
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        A0[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft0_term0
        A2[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/\
                (2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*\
                (sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

np.savetxt('data/A0_024w093.txt',sum_A0)
np.savetxt('data/A2_024w093.txt',sum_A2)
np.savetxt('data/G__024w093.txt',G)

