#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_93.txt') # LDS in perpendicular direction
#eiz_x = np.loadtxt('data/eiz_x_91.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in perpendicular direction
eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS of water, intervening medium
eiz_w = np.ones(len(eiz_x))#Aiz(eiz_x,eiz_z,eiz_w)
#eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w
#eiz_w = 1.001*eiz_x#np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#
#eiz_x = np.loadtxt('data/eiz_x_90.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_90.txt') # LDS in parallel direction
#
#eiz_x = np.loadtxt('data/eiz_x_91.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_91.txt') # LDS in parallel direction
#
#eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
#
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

Ls = np.arange(1e-9,10e-9,2e-9)  # separation distance between 2 cyclinders

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

a =  Aiz(eiz_x,eiz_z,eiz_w)
#a =  np.zeros(len(eiz_w))#Aiz(eiz_x,eiz_z,eiz_w)
#a =  2.* np.ones(len(eiz_w))#Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

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
        #print a,p,delta
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
                + 4.*(1.+a[j])*(1.+a[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
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
ratio =  sum_A0/sum_A2
percent= (sum_A0-sum_A2)/sum_A0
print 'para   = eiz65_z'
print 'perp   = H20'
print 'medium = H2O'
print 'A0 = ',sum_A0
print 'A2 = ',sum_A2
print 'A0/A2  = ',ratio
print '(A0-A2)/A0 = ',percent
#    print 1e9*Ls[i],1e21*(kbT/32)*A0[i,0]
#    print 1e9*Ls[i],1e21*(kbT/32)*A2[i,0]
#    print '-------------------'
#print sum_A0[0] 
#print sum_A2[0]
#print '******************'
#np.savetxt('data/A0_n_93.txt',A0)
#np.savetxt('data/A2_n_93.txt',A2)
#np.savetxt('data/Ls_n_93.txt',Ls)
#np.savetxt(  'data/A0_93_sum.txt',sum_A0)
#np.savetxt(  'data/A2_93_sum.txt',sum_A2)


#np.savetxt('data/A0_water_n_93.txt',A0)
#np.savetxt('data/A2_water_n_93.txt',A2)
#np.savetxt('data/A0_water_93_sum.txt',sum_A0)
#np.savetxt('data/A2_water_93_sum.txt',sum_A2)

#np.savetxt('data/A0_nonscreen_n_93.txt',A0)
#np.savetxt('data/A2_nonscreen_n_93.txt',A2)
#np.savetxt(  'data/A0_nonscreen_93_sum.txt',sum_A0)
#np.savetxt(  'data/A2_nonscreen_93_sum.txt',sum_A2)

#np.savetxt('data/A0_screen_n_93.txt',A0)
#np.savetxt('data/A2_screen_n_93.txt',A2)
#np.savetxt(  'data/A0_screen_93_sum.txt',sum_A0)
#np.savetxt(  'data/A2_screen_93_sum.txt',sum_A2)

#np.savetxt('data/A0_n0_equal_0_93.txt',A0)
#np.savetxt('data/A2_n0_equal_0_93.txt',A2)
#np.savetxt('data/Ls_n0_equal_0_93.txt',Ls)
#np.savetxt(  'data/A0_n0_equal_0_93_sum.txt',sum_A0)
#np.savetxt(  'data/A2_n0_equal_0_93_sum.txt',sum_A2)

#pl.figure()
#pl.plot(Ls,(kbT/32)*A0[:,  0], 'b:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  0]))
#pl.plot(Ls,(kbT/32)*A0[:,  1], 'g:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  1]))
#pl.plot(Ls,(kbT/32)*A0[:,  2], 'r:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  2]))
#pl.plot(Ls,(kbT/32)*A0[:,  3], 'c:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  3]))
#pl.plot(Ls,(kbT/32)*A0[:,  4], 'y:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  4]))
#pl.plot(Ls,(kbT/32)*A0[:,  5], 'm:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  5]))
#pl.plot(Ls,(kbT/32)*A0[:, 20], 'b-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 20]))
#pl.plot(Ls,(kbT/32)*A0[:, 40], 'g-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 40]))
#pl.plot(Ls,(kbT/32)*A0[:, 60], 'r-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 60]))
#pl.plot(Ls,(kbT/32)*A0[:, 80], 'c-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 80]))
#pl.plot(Ls,(kbT/32)*A0[:,100], 'b-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[100]))
#pl.plot(Ls,(kbT/32)*A0[:,200], 'g-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[200]))
#pl.plot(Ls,(kbT/32)*A0[:,300], 'r-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[300]))
#pl.plot(Ls,(kbT/32)*A0[:,400], 'c-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[400]))
##pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
#pl.legend(loc = 'best')      
##pl.title(r'93w93 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{n}\,\, [zJ]$')
#pl.xlabel(r'$\ell$')
#pl.savefig('plots/93_An0to400_vs_l.png')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()                    
#
#pl.figure()
#pl.plot(Ls,p[:,  0], 'b:' , label = r'$p(n=%2.1f)$'%(ns[  0]))
#pl.plot(Ls,p[:,  1], 'g:' , label = r'$p(n=%2.1f)$'%(ns[  1]))
#pl.plot(Ls,p[:,  2], 'r:' , label = r'$p(n=%2.1f)$'%(ns[  2]))
#pl.plot(Ls,p[:,  3], 'c:' , label = r'$p(n=%2.1f)$'%(ns[  3]))
#pl.plot(Ls,p[:,  4], 'y:' , label = r'$p(n=%2.1f)$'%(ns[  4]))
#pl.plot(Ls,p[:,  5], 'm:' , label = r'$p(n=%2.1f)$'%(ns[  5]))
#pl.plot(Ls,p[:, 20], 'b-.', label = r'$p(n=%2.1f)$'%(ns[ 20]))
#pl.plot(Ls,p[:, 40], 'g-.', label = r'$p(n=%2.1f)$'%(ns[ 40]))
#pl.plot(Ls,p[:, 60], 'r-.', label = r'$p(n=%2.1f)$'%(ns[ 60]))
#pl.plot(Ls,p[:, 80], 'c-.', label = r'$p(n=%2.1f)$'%(ns[ 80]))
#pl.plot(Ls,p[:,100], 'b-' , label = r'$p(n=%2.1f)$'%(ns[100]))
#pl.plot(Ls,p[:,200], 'g-' , label = r'$p(n=%2.1f)$'%(ns[200]))
#pl.plot(Ls,p[:,300], 'r-' , label = r'$p(n=%2.1f)$'%(ns[300]))
#pl.plot(Ls,p[:,400], 'c-' , label = r'$p(n=%2.1f)$'%(ns[400]))
##pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
#pl.legend(loc = 'best')      
##pl.title(r'65w65 Matsubara terms')
#pl.ylabel(r'$p_{n}(\ell)\,\, [zJ]$')
#pl.xlabel(r'$\ell$')
#pl.savefig('plots/93_pn0to400_vs_l.png')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()                    
#
#pl.figure()
#pl.plot(ns,p[ 0,:], 'b:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 0]))
#pl.plot(ns,p[ 1,:], 'g:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 1]))
#pl.plot(ns,p[ 2,:], 'r:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 2]))
#pl.plot(ns,p[ 3,:], 'c:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 3]))
#pl.plot(ns,p[ 4,:], 'y:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 4]))
#pl.plot(ns,p[ 5,:], 'm:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 5]))
#pl.plot(ns,p[10,:], 'b-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))
#pl.plot(ns,p[15,:], 'g-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[15]))
#pl.plot(ns,p[20,:], 'r-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
#pl.plot(ns,p[25,:], 'c-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[25]))
#pl.plot(ns,p[30,:], 'b-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
#pl.plot(ns,p[35,:], 'g-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[35]))
#pl.plot(ns,p[40,:], 'r-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
#pl.plot(ns,p[45,:], 'c-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[45]))
#pl.plot(ns, p[10,:]/p[10,:], 'k-')# , label = r'$p(l=%2.1f)$'%(Ls[45]))
#pl.legend(loc = 'best')      
##pl.title(r'65w65 p(n,l)')
#pl.ylabel(r'$p_{n}(\ell)$')
#pl.xlabel(r'$n$')
#pl.savefig('plots/93_pl_vs_n.pdf')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()                    

##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)

#        #A0[A0>1e6]= np.nan #NOTE: remove me later
#        #A2[A2>1e6]= np.nan #NOTE: remove me later
#    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
#    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
#pl.figure()
#pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
##pl.title(r'290w290 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
#pl.xlabel(r'$N$')
#pl.savefig('plots/65_A_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()
#
#pl.figure()
#pl.loglog(p[0,:],(kbT/(32.)) * A0[0,:] / A0[0,0], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(p[1,:],(kbT/(32.)) * A0[1,:] / A0[1,0], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(p[2,:],(kbT/(32.)) * A0[2,:] / A0[2,0], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(p[3,:],(kbT/(32.)) * A0[3,:] / A0[3,0], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.loglog(p[0,:],(kbT/(32.)) * A2[0,:] / A2[0,0], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(p[1,:],(kbT/(32.)) * A2[1,:] / A2[1,0], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(p[2,:],(kbT/(32.)) * A2[2,:] / A2[2,0], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(p[3,:],(kbT/(32.)) * A2[3,:] / A2[3,0], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms vs p')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
##pl.title(r'290w290 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
#pl.xlabel(r'$N$')
#pl.savefig('plots/65_A_vs_p.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()
#
#print 'A0(separation) = ',sum_A0
#print 'A2(separation) = ',sum_A2
#print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
#print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
#
#np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
#
A_py_par  = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
A0_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
A2_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
A_GH_par  = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'

x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
y_ax_per = r'$\mathrm{\mathcal{A}_\perp (\ell)}\,\,\,\rm{[zJ]}$'

def title(cnt1,cnt2,orientation):
	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)

def svfig(cnt1,cnt2,orientation):
	return 'plots/140322_%sw%s_HCs_%s_ret.png'%(cnt1,cnt2,orientation)

pl.figure()
pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.legend(loc = 'best')
#
#pl.title(title('6','5','Log-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
#pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
#pl.savefig(svfig('65','65','perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
pl.savefig(svfig('93','93','n0_equal_0_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
pl.show()
#
pl.figure()
pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'%(1e9*Ls[0],      1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('6','5','perpendicular'))
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.legend(loc = 'best')
#pl.title(title('6','5','Semi-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
#pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
#pl.savefig(svfig('65','65','semilog_perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
pl.savefig(svfig('93','93','n0_equal_0_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
pl.show()
#
#





