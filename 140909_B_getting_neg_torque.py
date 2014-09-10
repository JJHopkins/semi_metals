#!/usr/bin/env python
import sys
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import trapz,romb,quad
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show
from pylab import *
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

r_1 = 0.5e-9 # ??? XXX DNA?
#r_1 = 0.4234e-9 # 93
r_2 = 0.4234e-9

filename0=sys.argv[0] # eiz_x
print filename0

eiz_x_1 = np.loadtxt('data/eiz_x_CG.txt') # LDS in perpendicular direction
eiz_z_1 = np.loadtxt('data/eiz_z_CG.txt') # LDS in parallel direction

Eiz_x_2 = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
Eiz_z_2 = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction

eiz_x_2 =Eiz_x_2[:500]   
eiz_z_2 =Eiz_z_2[:500]  


eiz_w = np.loadtxt('data/eiz_i_Al.txt') # LDS of water, intervening medium
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium
medium = 'Al2O3'

# Constants
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]
# Matsubara frequencies
ns = np.arange(0.,500) 
zs = ns * coeff         
#ls = [100e-9]
#Ls = np.arange(1e-9,151e-9,5e-9)  # separation distance between 2 cyclinders
R1 = 0.5e-9
Ls = np.linspace(1.1*R1,5e-9,10)
#Ls = np.linspace(1e-9,199e-9,41)
#ls = np.arange(1e-9,151e-9,5e-9)  # separation distance between 2 cyclinders

#thetas  = [np.pi/8,0.780,0.785,0.790]#np.linspace(0.01,np.pi-0.01,17)
thetas  = np.linspace(0.01,np.pi/2-0.01,9)

#thetas= [np.pi/8,np.pi/4]#np.linspace(0.775,0.795,3)
#thetas  = np.linspace(0.775,0.795,3)

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

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))
Tau = np.zeros(shape = (len(Ls),len(thetas)))

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta_par = Delta_par(eiz_z,eiz_w)
delta_prp = Delta_prp(eiz_x,eiz_w)
 
# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a[0])*(1.+3.*a[0])
        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a[0])*(1.-a[0])
        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
Ft0_term0 = romb(f0_term0)
#Ft0_term0 = np.sum(f0_term0)
Ft2_term0 = romb(f2_term0)
#Ft2_term0 = np.sum(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *2.*((delta_par[j]+6.*delta_prp[j])*(delta_par[j]+6.*delta_prp[j])*T*T*T*T\
                +2.*(delta_par[j]*delta_par[j]+4.*delta_prp[j]*delta_par[j]+4.*delta_prp[j]*delta_par[j]+12.*delta_prp[j]*delta_prp[j])*T*T \
                + 2.*(delta_par[j]+2.*delta_prp[j])*(delta_par[j]+2.*delta_prp[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *(T*T*T*T+4.*(T*T)+4.)*(delta_par[j]-2.*delta_prp[j])*(delta_par[j]-2.*delta_prp[j])
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        Ft0 = romb(f0)
        Ft2 = romb(f2)

#        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
#                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
#                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
#                + 4.*(1.+a[j])*(1.+a[j]))
#                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
#        # Integrand A2                
#        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
#                /(T*T+1.)\
#                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))
#                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
#        Ft0 = romb(f0)
#        Ft2 = romb(f2)

        A0[i,j] = p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,j] = delta_par[j]*delta_par[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,0] = (1./2) * delta_prp[0]*delta_prp[0]*Ft0_term0
        A0[i,0] = 0.#(1./2) * delta_prp[0]*delta_prp[0]*Ft0_term0
        
        A2[i,j] = p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,j] = delta_par[j]*delta_par[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,0] = (1./2) * delta_prp[0]*delta_prp[0]*Ft2_term0
        A2[i,0] = 0.#(1./2) * delta_prp[0]*delta_prp[0]*Ft2_term0
    #sum_A0 = kbT * np.sum(A0, axis = 1)
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    #sum_A2 = kbT * np.sum(A2, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for j,theta in enumerate(thetas):
        G[i,j]=-(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[j])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[j]))
        Tau[i,j]=-(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*Ls[i]*Ls[i]*Ls[i]*Ls[i])\
                *(sum_A0[i]+sum_A2[i]*(2.-np.cos(2.*thetas[j])))\
                *np.cos(thetas[j])/(np.sin(thetas[j])*np.sin(thetas[j]))
#np.savetxt('data/A0_cyl_93w93.txt',sum_A0)
#np.savetxt('data/A2_cyl_93w93.txt',sum_A2)
#np.savetxt('data/G_cyl_vs_R2_93w93_screened.txt',G)
#np.savetxt('data/T_cyl_93w93.txt',Tau)
#sys.exit()
pl.figure()
pl.loglog(1e9*Ls,-1e21*Tau[:,2],'b-',label = r'$M(\ell,\pi/8)_{cyl}$')
#pl.loglog(1e9*ls, 1e21*T[:,2], 'g-',label = r'$M(\ell,\pi/8)_{plane}$')
pl.title(r'Torque at $\pi$/8 vs separation [9,3]w[9,3]')
#pl.xlabel(r'$angle\,\,\,[radians]$')
pl.xlabel(r'$separation\,\,\,[nm]$')
pl.ylabel(r'$M(\ell,\theta)\,\,\,[zJ]$')
#pl.axis([0,180,0.3,3.0])
pl.legend(loc = 'best')
#pl.savefig('plots/loglog_Tcp_v_l.png')
pl.show()

pl.plot(ns, eiz_w    ,'k:',label=r'$\epsilon_{h20}$')
pl.plot(ns, eiz_z    ,'b:',label=r'$\epsilon_{\parallel}$')
pl.plot(ns,delta_par, 'b-',label=r'$\Delta_{\parallel}\,max=%2.2f$'%(np.max(delta_par)))
pl.plot(ns, eiz_x    ,'r:',label=r'$\epsilon_{\perp}$')
pl.plot(ns,delta_prp ,'r-',label=r'$\Delta_{\perp}\,max=%2.2f$'%(np.max(delta_prp)))
pl.plot(ns,2.*delta_prp ,'r--',label=r'$2\Delta_{\perp}$')
pl.plot(ns,delta_par-2.*delta_prp,'c',label=r'$\Delta_{\parallel}-2\Delta_{\perp}\,max=%2.2f$'%(np.max(delta_par-2.*delta_prp)))
pl.plot(ns,(delta_par-2.*delta_prp)**2.,'m',label=r'$(\Delta_{\parallel}-2\Delta_{\perp})^2\,max=%2.2f$'%(np.max((delta_par-2.*delta_prp))**2))
#pl.plot(ns, a,'y', label=r'$a(i\omega)$')
pl.legend(loc = 'best')
pl.xlabel(r'$n$')
#pl.ylabel(r'$n$')
#pl.savefig('plots/eiz_delta_a_check.pdf')
pl.show()
sys.exit()
#
#pl.figure()
#pl.loglog(ns, eiz_w    ,'k:',linewidth = 2,label=r'$\epsilon_{h20}$')
#pl.loglog(ns, eiz_z    ,'b:',linewidth = 2,label=r'$\epsilon_{\parallel}$')
#pl.loglog(ns,delta_par, 'b-',linewidth = 2,label=r'$\Delta_{\parallel}\,max=%2.2f$'%(np.max(delta_par)))
#pl.loglog(ns, eiz_x    ,'r:',linewidth = 2,label=r'$\epsilon_{\perp}$')
#pl.loglog(ns,delta_prp ,'r-',linewidth = 2,label=r'$\Delta_{\perp}\,max=%2.2f$'%(np.max(delta_prp)))
#pl.loglog(ns,2.*delta_prp ,'r--',linewidth = 2,label=r'$2\Delta_{\perp}$')
#pl.loglog(ns,np.abs(delta_par-2.*delta_prp),'c',linewidth=2,label=r'$|\Delta_{\parallel}-2\Delta_{\perp}|\,\,max=%2.2f$'%(np.max(delta_par-2.*delta_prp)))
#pl.loglog(ns,(delta_par-2.*delta_prp)**2.,'m',linewidth = 2,label=r'$(\Delta_{\parallel}-2\Delta_{\perp})^2\,max=%2.2f$'%(np.max((delta_par-2.*delta_prp))**2))
##pl.plot(ns, a,'y', label=r'$a(i\omega)$')
#pl.legend(loc = 'best')
#pl.xlabel(r'$n$')
##pl.ylabel(r'$n$')
#pl.savefig('plots/loglog_eiz_delta_a_check.pdf')
#pl.show()
#
#pl.figure()
#pl.semilogx(ns, eiz_w    ,'k:',linewidth = 2,label=r'$\epsilon_{h20}$')
#pl.semilogx(ns, eiz_z    ,'b:',linewidth = 2,label=r'$\epsilon_{\parallel}$')
#pl.semilogx(ns,delta_par, 'b-',linewidth = 2,label=r'$\Delta_{\parallel}\,max=%2.2f$'%(np.max(delta_par)))
#pl.semilogx(ns, eiz_x    ,'r:',linewidth = 2,label=r'$\epsilon_{\perp}$')
#pl.semilogx(ns,delta_prp ,'r-',linewidth = 2,label=r'$\Delta_{\perp}\,max=%2.2f$'%(np.max(delta_prp)))
#pl.semilogx(ns,2.*delta_prp ,'r--',linewidth = 2,label=r'$2\Delta_{\perp}$')
#pl.semilogx(ns,delta_par-2.*delta_prp,'c',linewidth = 2,label=r'$\Delta_{\parallel}-2\Delta_{\perp}\,max=%2.2f$'%(np.max(delta_par-2.*delta_prp)))
#pl.semilogx(ns,(delta_par-2.*delta_prp)**2.,'m',linewidth = 2,label=r'$(\Delta_{\parallel}-2\Delta_{\perp})^2\,max=%2.2f$'%(np.max((delta_par-2.*delta_prp))**2))
##pl.plot(ns, a,'y', label=r'$a(i\omega)$')
#pl.legend(loc = 'best')
#pl.xlabel(r'$n$')
##pl.ylabel(r'$n$')
#pl.savefig('plots/semilogx_eiz_delta_a_check.pdf')
#pl.show()
#
##pl.figure()
##pl.plot(ns, delta_par,'b:')#, label=r'$screened$')
##pl.plot(ns, delta_prp,'g:')#, label=r'$screened \Delta_{\perp}$')
##pl.plot(ns, delta_prp,'r:')#, label=r'$screened \Delta_{\parallel}-2\Delta_{\perp}$')
##pl.plot(ns, a,'c:')#, label=r'$screened a(i\omega)$')
##pl.legend(loc = 'best')
##pl.xlabel(r'$n$')
###pl.ylabel(r'$n$')
##pl.show()
#
##    print 1e9*Ls[i],1e21*(kbT/32)*A0[i,0]
##    print 1e9*Ls[i],1e21*(kbT/32)*A2[i,0]
##    print '-------------------'
##print sum_A0[0] 
##print sum_A2[0]
##print '******************'
##np.savetxt('data/A0_n_93.txt',A0)
##np.savetxt('data/A2_n_93.txt',A2)
##np.savetxt('data/Ls_n_93.txt',Ls)
##np.savetxt(  'data/A0_93_sum.txt',sum_A0)
##np.savetxt(  'data/A2_93_sum.txt',sum_A2)
#
#np.savetxt('data/A0_cyl_93w93.txt',sum_A0)
#np.savetxt('data/A2_cyl_93w93.txt',sum_A2)
#
#
##np.savetxt('data/A0_water_n_93.txt',A0)
##np.savetxt('data/A2_water_n_93.txt',A2)
##np.savetxt('data/A0_water_93_sum.txt',sum_A0)
##np.savetxt('data/A2_water_93_sum.txt',sum_A2)
#
##np.savetxt('data/A0_nonscreen_n_93.txt',A0)
##np.savetxt('data/A2_nonscreen_n_93.txt',A2)
##np.savetxt(  'data/A0_nonscreen_93_sum.txt',sum_A0)
##np.savetxt(  'data/A2_nonscreen_93_sum.txt',sum_A2)
#
##np.savetxt('data/A0_screen_n_93.txt',A0)
##np.savetxt('data/A2_screen_n_93.txt',A2)
##np.savetxt(  'data/A0_screen_93_sum.txt',sum_A0)
##np.savetxt(  'data/A2_screen_93_sum.txt',sum_A2)
#
##np.savetxt('data/A0_n0_equal_0_93.txt',A0)
##np.savetxt('data/A2_n0_equal_0_93.txt',A2)
##np.savetxt('data/Ls_n0_equal_0_93.txt',Ls)
##np.savetxt(  'data/A0_n0_equal_0_93_sum.txt',sum_A0)
##np.savetxt(  'data/A2_n0_equal_0_93_sum.txt',sum_A2)
#
##pl.figure()
##pl.plot(Ls,(kbT/32)*A0[:,  0], 'b:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  0]))
##pl.plot(Ls,(kbT/32)*A0[:,  1], 'g:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  1]))
##pl.plot(Ls,(kbT/32)*A0[:,  2], 'r:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  2]))
##pl.plot(Ls,(kbT/32)*A0[:,  3], 'c:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  3]))
##pl.plot(Ls,(kbT/32)*A0[:,  4], 'y:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  4]))
##pl.plot(Ls,(kbT/32)*A0[:,  5], 'm:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  5]))
##pl.plot(Ls,(kbT/32)*A0[:, 20], 'b-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 20]))
##pl.plot(Ls,(kbT/32)*A0[:, 40], 'g-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 40]))
##pl.plot(Ls,(kbT/32)*A0[:, 60], 'r-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 60]))
##pl.plot(Ls,(kbT/32)*A0[:, 80], 'c-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 80]))
##pl.plot(Ls,(kbT/32)*A0[:,100], 'b-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[100]))
##pl.plot(Ls,(kbT/32)*A0[:,200], 'g-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[200]))
##pl.plot(Ls,(kbT/32)*A0[:,300], 'r-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[300]))
##pl.plot(Ls,(kbT/32)*A0[:,400], 'c-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[400]))
###pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
##pl.legend(loc = 'best')      
###pl.title(r'93w93 Matsubara terms')
##pl.ylabel(r'$\mathcal{A}^{(0)}_{n}\,\, [zJ]$')
##pl.xlabel(r'$\ell$')
##pl.savefig('plots/93_An0to400_vs_l.png')
####pl.savefig('plots/90_A_vs_n.pdf')
####pl.savefig('plots/91_A_vs_n.pdf')
####pl.savefig('plots/93_A_vs_n.pdf')
####pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()                    
##
##pl.figure()
##pl.plot(Ls,p[:,  0], 'b:' , label = r'$p(n=%2.1f)$'%(ns[  0]))
##pl.plot(Ls,p[:,  1], 'g:' , label = r'$p(n=%2.1f)$'%(ns[  1]))
##pl.plot(Ls,p[:,  2], 'r:' , label = r'$p(n=%2.1f)$'%(ns[  2]))
##pl.plot(Ls,p[:,  3], 'c:' , label = r'$p(n=%2.1f)$'%(ns[  3]))
##pl.plot(Ls,p[:,  4], 'y:' , label = r'$p(n=%2.1f)$'%(ns[  4]))
##pl.plot(Ls,p[:,  5], 'm:' , label = r'$p(n=%2.1f)$'%(ns[  5]))
##pl.plot(Ls,p[:, 20], 'b-.', label = r'$p(n=%2.1f)$'%(ns[ 20]))
##pl.plot(Ls,p[:, 40], 'g-.', label = r'$p(n=%2.1f)$'%(ns[ 40]))
##pl.plot(Ls,p[:, 60], 'r-.', label = r'$p(n=%2.1f)$'%(ns[ 60]))
##pl.plot(Ls,p[:, 80], 'c-.', label = r'$p(n=%2.1f)$'%(ns[ 80]))
##pl.plot(Ls,p[:,100], 'b-' , label = r'$p(n=%2.1f)$'%(ns[100]))
##pl.plot(Ls,p[:,200], 'g-' , label = r'$p(n=%2.1f)$'%(ns[200]))
##pl.plot(Ls,p[:,300], 'r-' , label = r'$p(n=%2.1f)$'%(ns[300]))
##pl.plot(Ls,p[:,400], 'c-' , label = r'$p(n=%2.1f)$'%(ns[400]))
###pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
##pl.legend(loc = 'best')      
###pl.title(r'65w65 Matsubara terms')
##pl.ylabel(r'$p_{n}(\ell)\,\, [zJ]$')
##pl.xlabel(r'$\ell$')
##pl.savefig('plots/93_pn0to400_vs_l.png')
####pl.savefig('plots/90_A_vs_n.pdf')
####pl.savefig('plots/91_A_vs_n.pdf')
####pl.savefig('plots/93_A_vs_n.pdf')
####pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()                    
##
##pl.figure()
##pl.plot(ns,p[ 0,:], 'b:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 0]))
##pl.plot(ns,p[ 1,:], 'g:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 1]))
##pl.plot(ns,p[ 2,:], 'r:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 2]))
##pl.plot(ns,p[ 3,:], 'c:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 3]))
##pl.plot(ns,p[ 4,:], 'y:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 4]))
##pl.plot(ns,p[ 5,:], 'm:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 5]))
##pl.plot(ns,p[10,:], 'b-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))
##pl.plot(ns,p[15,:], 'g-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[15]))
##pl.plot(ns,p[20,:], 'r-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
##pl.plot(ns,p[25,:], 'c-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[25]))
##pl.plot(ns,p[30,:], 'b-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
##pl.plot(ns,p[35,:], 'g-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[35]))
##pl.plot(ns,p[40,:], 'r-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
##pl.plot(ns,p[45,:], 'c-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[45]))
##pl.plot(ns, p[10,:]/p[10,:], 'k-')# , label = r'$p(l=%2.1f)$'%(Ls[45]))
##pl.legend(loc = 'best')      
###pl.title(r'65w65 p(n,l)')
##pl.ylabel(r'$p_{n}(\ell)$')
##pl.xlabel(r'$n$')
##pl.savefig('plots/93_pl_vs_n.pdf')
####pl.savefig('plots/90_A_vs_n.pdf')
####pl.savefig('plots/91_A_vs_n.pdf')
####pl.savefig('plots/93_A_vs_n.pdf')
####pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()                    
#
###np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
#
##        #A0[A0>1e6]= np.nan #NOTE: remove me later
##        #A2[A2>1e6]= np.nan #NOTE: remove me later
##    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
##    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
##pl.figure()
##pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.legend(loc = 'best')
##pl.title(r'65w65 Matsubara terms')
###pl.title(r'90w90 Matsubara terms')
###pl.title(r'91w91 Matsubara terms')
###pl.title(r'93w93 Matsubara terms')
###pl.title(r'290w290 Matsubara terms')
##pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
##pl.xlabel(r'$N$')
##pl.savefig('plots/65_A_vs_n.pdf')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()
##
##pl.figure()
##pl.loglog(p[0,:],(kbT/(32.)) * A0[0,:] / A0[0,0], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(p[1,:],(kbT/(32.)) * A0[1,:] / A0[1,0], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(p[2,:],(kbT/(32.)) * A0[2,:] / A0[2,0], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(p[3,:],(kbT/(32.)) * A0[3,:] / A0[3,0], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.loglog(p[0,:],(kbT/(32.)) * A2[0,:] / A2[0,0], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(p[1,:],(kbT/(32.)) * A2[1,:] / A2[1,0], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(p[2,:],(kbT/(32.)) * A2[2,:] / A2[2,0], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(p[3,:],(kbT/(32.)) * A2[3,:] / A2[3,0], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.legend(loc = 'best')
##pl.title(r'65w65 Matsubara terms vs p')
###pl.title(r'90w90 Matsubara terms')
###pl.title(r'91w91 Matsubara terms')
###pl.title(r'93w93 Matsubara terms')
###pl.title(r'290w290 Matsubara terms')
##pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
##pl.xlabel(r'$N$')
##pl.savefig('plots/65_A_vs_p.pdf')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()
##
##print 'A0(separation) = ',sum_A0
##print 'A2(separation) = ',sum_A2
##print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
##print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
##
##np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
##
#A_py_par  = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
#A0_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
#A2_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
#A_GH_par  = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'
#
#x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
#y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
#y_ax_per = r'$\mathrm{\mathcal{A}_\perp (\ell)}\,\,\,\rm{[zJ]}$'
#
#def title(cnt1,cnt2,orientation):
#	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)
#
#def svfig(cnt1,cnt2,orientation):
#	return 'plots/140322_%sw%s_HCs_%s_ret.png'%(cnt1,cnt2,orientation)
#
#pl.figure()
#pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
#pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.legend(loc = 'best')
##
##pl.title(title('6','5','Log-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
##pl.savefig(svfig('65','65','perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','n0_equal_0_perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
##
#pl.figure()
#pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'%(1e9*Ls[0],      1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
#pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.title(title('6','5','perpendicular'))
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.legend(loc = 'best')
##pl.title(title('6','5','Semi-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
##pl.savefig(svfig('65','65','semilog_perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','n0_equal_0_perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
##
##
#
#
#
#
#


