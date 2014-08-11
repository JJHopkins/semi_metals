#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

filename0=sys.argv[0]
print filename0

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
#eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
#eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#eiz_x = eiz_X[:500]
#eiz_z = eiz_Z[:500]
eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
#eiz_w  = np.ones(len(eiz_x))

r_1 = 1.6272e-9 #24
r_2 = 1.6272e-9
#r_1 = 0.3734e-9 #65
#r_2 = 0.3734e-9
#
#r_1 = 0.4234e-9 # 93
#r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))
NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
#print '93v93, unscreened'
#print '93w93, screened'
#print '65v65, unscreened'
#print '65w65, screened'
print '24w24, screened'
#print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.figure()
pl.loglog(Ls,1e21*G[:,2] ,'b-',label= '24w24')
pl.loglog(Ls,-1e12*NN[:,2],'b:')


#------------------------------------------------------------------------------
# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
#eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
#eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#eiz_x = eiz_X[:500]
#eiz_z = eiz_Z[:500]
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
eiz_w  = np.ones(len(eiz_x))

r_1 = 1.6272e-9 #24
r_2 = 1.6272e-9
#r_1 = 0.3734e-9 #65
#r_2 = 0.3734e-9
#
#r_1 = 0.4234e-9 # 93
#r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
#print '93v93, unscreened'
#print '93w93, screened'
#print '65v65, unscreened'
#print '65w65, screened'
#print '24w24, screened'
print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.loglog(Ls,1e21*G[:,2] ,'g-',label= '24v24')
pl.loglog(Ls,-1e12*NN[:,2],'g:')


#65w65 screen------------------------------------------------------------------------------
# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
#eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#eiz_x = eiz_X[:500]
#eiz_z = eiz_Z[:500]
eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
#eiz_w  = np.ones(len(eiz_x))

#r_1 = 1.6272e-9 #24
#r_2 = 1.6272e-9
r_1 = 0.3734e-9 #65
r_2 = 0.3734e-9
#
#r_1 = 0.4234e-9 # 93
#r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
#print '93v93, unscreened'
#print '93w93, screened'
#print '65v65, unscreened'
print '65w65, screened'
#print '24w24, screened'
#print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.loglog(Ls,1e21*G[:,2] ,'r-',label= '65w65')
pl.loglog(Ls,-1e12*NN[:,2],'r:')

#------------------------------------------------------------------------------
# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
#eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#eiz_x = eiz_X[:500]
#eiz_z = eiz_Z[:500]
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
eiz_w  = np.ones(len(eiz_x))

#r_1 = 1.6272e-9 #24
#r_2 = 1.6272e-9
r_1 = 0.3734e-9 #65
r_2 = 0.3734e-9
#
#r_1 = 0.4234e-9 # 93
#r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
#print '93v93, unscreened'
#print '93w93, screened'
print '65v65, unscreened'
#print '65w65, screened'
#print '24w24, screened'
#print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.loglog(Ls,1e21*G[:,2] ,'c-',label= '65v65')
pl.loglog(Ls,-1e12*NN[:,2],'c:')



#93w93 screen------------------------------------------------------------------------------
# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
#eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
eiz_x = eiz_X[:500]
eiz_z = eiz_Z[:500]
eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
#eiz_w  = np.ones(len(eiz_x))

#r_1 = 1.6272e-9 #24
#r_2 = 1.6272e-9
#r_1 = 0.3734e-9 #65
#r_2 = 0.3734e-9
#
r_1 = 0.4234e-9 # 93
r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
#print '93v93, unscreened'
print '93w93, screened'
#print '65v65, unscreened'
#print '65w65, screened'
#print '24w24, screened'
#print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.loglog(Ls,1e21*G[:,2] ,'m-',label= '93w93')
pl.loglog(Ls,-1e12*NN[:,2],'m:')


#------------------------------------------------------------------------------
# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
#eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
eiz_X = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
eiz_Z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
eiz_x = eiz_X[:500]
eiz_z = eiz_Z[:500]
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
eiz_w  = np.ones(len(eiz_x))

#r_1 = 1.6272e-9 #24
#r_2 = 1.6272e-9
#r_1 = 0.3734e-9 #65
#r_2 = 0.3734e-9
#
r_1 = 0.4234e-9 # 93
r_2 = 0.4234e-9

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(10e-9,40e-9,1e-9)  #separation between 2 cyclinders
thetas = np.arange(0.783,0.786,0.001)

#Integration vars
T  = np.linspace(0.,2.**17, 1.+2.**17)
U  = np.linspace(0.,2.**17, 1.+2.**17)

# Define functions
def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])\
        *((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    #print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))\
                *((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))
        #Ft0 = np.sum(f0)
        #Ft2 = np.sum(f2)
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Ft = romb(f , axis = 1)
        #Fty =romb(Ft, axis = 0)

        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))

NN,TT = np.gradient(G,Ls[1]-Ls[0],thetas[1]-thetas[0])
#T = -np.diff(G[2,:])/np.diff(thetas)
#N = -np.diff(G[:,2])/np.diff(Ls)
print'   '
print'-----------------------'
print '93v93, unscreened'
#print '93w93, screened'
#print '65v65, unscreened'
#print '65w65, screened'
#print '24w24, screened'
#print '24v24, unscreened'
print'-----------------------'
print 'theta =', thetas[2]
print 'Ls    = ',Ls[2]
print 'A0(5nm) = ',1e21*sum_A0[2] 
print 'A2(5nm) = ',1e21*sum_A2[2]
print 'G  = ',1e21*G[2,2]
print 'T  = ',-1e21*TT[2,2]
print 'N  = ',-1e12*NN[2,2]

pl.loglog(Ls,1e21*G[:,2] ,'y-',label= '93v93')
pl.loglog(Ls,-1e12*NN[:,2],'y:')
pl.legend(loc = 'best')
pl.title(r'Free Energy and Normal Force')
pl.xlabel('separation [nm]')
pl.ylabel('$G\,\,[zJ],\,\,\,\,|F_{norm}|\,\,[pN]$')
pl.savefig('plots/G_N_for_DD_test_cases.pdf')
pl.show()



