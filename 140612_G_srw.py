#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

#sum_A0 = np.loadtxt('data/A0_51w290_perpendicular_ret.txt')
#sum_A2 = np.loadtxt('data/A2_51w290_perpendicular_ret.txt')
sum_A0 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_65_sum.txt')
sum_A2 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_65_sum.txt')
#sum_A0 = np.loadtxt('data/A0_n0_equal_0_93_sum.txt')
#sum_A2 = np.loadtxt('data/A2_n0_equal_0_93_sum.txt')
#sum_A0 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A0_93_sum.txt')
#sum_A2 = np.loadtxt('data/Completed_data_n0to500_499_L1to450_449/A2_93_sum.txt')

#eiz_x_051 = np.loadtxt('data/eiz_x_51.txt') # LDS in perpendicular direction
#eiz_z_051 = np.loadtxt('data/eiz_z_51.txt') # LDS in parallel direction
#
#eiz_x_290 = np.loadtxt('data/eiz_x_290.txt') # LDS in perpendicular direction
#eiz_z_290 = np.loadtxt('data/eiz_z_290.txt') # LDS in parallel direction
#
#eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
##eiz_w[0] = 79.0

# Constants
#r_1 = 0.5e-9
r_1 = 0.4234e-9
#r_2 = 0.5e-9
r_2 = 0.432e-9
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
#coeff = 0.159           # [eV]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
#ns = z/(0.159)#

ns = np.arange(0.,500.) 
zs = ns * coeff         

#Ls = np.arange(1e-9,100e-9,1e-9)  # separation distance between 2 cyclinders
thetas = np.arange(np.pi/12, np.pi/2,0.001)
Ls = np.arange(1e-9,450e-9,1e-9)  # separation distance between 2 cyclinders

G = np.zeros(shape = (len(Ls),len(thetas)))

for i,L in enumerate(Ls):
    for j,theta in enumerate(thetas):
        G[i,j]=(-np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))
        #G[i,j]=(-np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*1e9*1e9*1e9*1e9*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(1e21*sum_A0[i]+1e21*sum_A2[i]*np.cos(2.*thetas[j]))
        #print i,j,G
#    pl.plot(ns,A2_org[i,:],'-')
#    pl.plot(ns,-A2_new[i,:],'--')
#pl.show()

dGdl  = np.diff(G, axis = 0)
dGdth = np.diff(G, axis = 1)
lForce =  np.log(dGdl)
lTau   = np.log(dGdth)
#print 1e21*dGdl[4,52]  
#print 1e23*dGdth[4,52] 
#
#print lForce[4,52]
#print lTau[4,52]
#Force = -dGdl
#Tau   = -dGdth

#pl.plot(1e9*Ls[:98],np.diff(1e21*sum_A0),      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[34],1e21*sum_A0[34]))
#pl.figure()
#pl.loglog(Ls[:-1],1e12*1e9*dGdl[:,10], label=r'$\theta =%1.2f rad$'%(thetas[10]))
#pl.loglog(Ls[:-1],1e12*1e9*dGdl[:,52], label=r'$\theta =%1.2f rad$'%(thetas[52]))
#pl.loglog(Ls[:-1],1e12*1e9*dGdl[:,104], label=r'$\theta=%1.2f rad$'%(thetas[104]))
#pl.legend()
#pl.xlabel('[meters]')
#pl.ylabel(r'$dG/d\ell\,\,[pN]$')
##pl.savefig('dGdlpN.pdf')
#pl.show()
#
pl.figure()
pl.loglog(Ls,-1e21*G[:,10], label=r'$\theta =%1.2f rad$'%(thetas[10]))
pl.loglog(Ls,-1e21*G[:,52], label=r'$\theta =%1.2f rad$'%(thetas[52]))
pl.loglog(Ls,-1e21*G[:,104], label=r'$\theta=%1.2f rad$'%(thetas[104]))
pl.legend()
pl.xlabel('[meters]')
pl.ylabel(r'$-G\,\,[J]$')
#pl.savefig('GpN.pdf')
pl.show()
#
#pl.figure()
#pl.loglog(thetas[:-1],1e24*dGdth[3,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[3]))
#pl.loglog(thetas[:-1],1e24*dGdth[4,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[4]))
#pl.loglog(thetas[:-1],1e24*dGdth[5,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[5]))
#pl.legend()
#pl.xlabel('[radians]')
#pl.ylabel(r'$dG/d\theta\,\,[J]/radian$')
##pl.savefig('dGdthpN.pdf')
#pl.show()
#
#pl.figure()
#pl.semilogy((180./np.pi)*thetas[:-1],1e24*dGdth[3,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[3]))
#pl.semilogy((180./np.pi)*thetas[:-1],1e24*dGdth[4,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[4]))
#pl.semilogy((180./np.pi)*thetas[:-1],1e24*dGdth[5,:], label=r'$\ell =%1.1f nm$'%(1e9*Ls[5]))
#pl.legend()
#pl.xlabel('[degrees]')
#pl.ylabel(r'$dG/d\theta\,\,[zJ]/angle$')
##pl.savefig('dGdthsemilogpN.pdf')
#pl.show()

#pl.figure()
#pl.plot(Ls[:-1],lForce[:,10], label=r'$\theta =%1.1f rad$'%(thetas[10]))
#pl.plot(Ls[:-1],lForce[:,52], label=r'$\theta =%1.1f rad$'%(thetas[52]))
#pl.plot(Ls[:-1],lForce[:,104], label=r'$\theta =%1.1f rad$'%(thetas[104]))
#pl.show()
#
#pl.figure()
#pl.plot(thetas[:-1],lTau[3,:])
#pl.plot(thetas[:-1],lTau[4,:])
#pl.plot(thetas[:-1],lTau[5,:])
#pl.show()
#
#
#pl.figure()
#pl.loglog(1e9*Ls,-1e21*G[:,1],'b--')
#pl.show()
#
#
#pl.figure()
#pl.loglog(1e9*Ls,-1e21*G[:,1],'b--')
#pl.show()
#
#pl.figure()
#pl.loglog(1e9*Ls,G[:,1]/G[0,1],'b--')
#pl.show()
#
