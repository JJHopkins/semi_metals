#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

filename0=sys.argv[0]
print filename0

eiz_x_1 = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z_1 = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction

#eiz_x_1 = np.loadtxt('data/eiz_x_CG.txt') # LDS in perpendicular direction
#eiz_z_1 = np.loadtxt('data/eiz_z_CG.txt') # LDS in parallel direction

Eiz_x_2 = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
Eiz_z_2 = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction

eiz_x_2 =Eiz_x_2[:500]   
eiz_z_2 =Eiz_z_2[:500]  

#eiz_x_1 =Eiz_x_2[:500]   
#eiz_z_1 =Eiz_z_2[:500]  

eiz_w = np.loadtxt('data/eiz_i_Al.txt') # LDS of water, intervening medium
#eiz_w = np.loadtxt('data/140730_eiz_w.txt') # LDS of water, intervening medium

# Constants
r_1 = 0.5e-9
r_2 = 0.5e-9
#thetas = [0.78,0.785,0.79]#np.pi/2
thetas = np.linspace(0.1,3.1,301)
#thetas = np.linspace(0.1,1.50,141)
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(1e-9,106e-9,5e-9)  # separation distance between 2 cyclinders

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

#def DDDg2(perp1, perp2, par1,par2, med):
#    return (2./med)*((perp1-med)*(perp2-med)/(2.*med) -((perp1+med)*(par1-med)*(perp2-med)-(perp2+med)*(par2-med)*(perp1-med)+2.*med*(perp1-med)*(perp2-med))/((perp2+med)*(perp1+med)))
#
#def DDg2(perp1, perp2, par1,par2, med):
#    return (1./(med*med))*(par2-med)*(par1-med)*\
#            ((par1-med)/med-2.*(perp1-med)/(perp1+med))*\
#            ((par2-med)/med -2.*(perp2-med)/(perp2+med))

def Pn(e,zn,l):
	return np.sqrt(e)*zn*l*(1./c)

p = np.zeros(shape = (len(Ls),len(ns)))
A0 = np.zeros(shape = (len(Ls),len(ns)))
A2 = np.zeros(shape = (len(Ls),len(ns)))
A2_org = np.zeros(shape = (len(Ls),len(ns)))
A2_new = np.zeros(shape = (len(Ls),len(ns)))
G = np.zeros(shape = (len(Ls),len(thetas)))
Tau = np.zeros(shape = (len(Ls),len(thetas)))

a_1 =   Aiz(eiz_x_1,eiz_z_1,eiz_w)
a_2 =   Aiz(eiz_x_2,eiz_z_2,eiz_w)
delta_1 = Delta(eiz_z_1,eiz_w)
delta_2 = Delta(eiz_z_2,eiz_w)
delta_prp1 = Delta_perp(eiz_x_1,eiz_w)
delta_prp2 = Delta_perp(eiz_x_2,eiz_w)
#D   = delta_prp2*delta_1*(delta_1 - 2.*delta_prp1)*(delta_2-2.*delta_prp2)
#DD  = DDDg2(eiz_x_1,eiz_x_2,eiz_z_1,eiz_z_2,eiz_w)
#DDD = DDDg2(eiz_x_1,eiz_x_2,eiz_z_1,eiz_z_2,eiz_w)

# Integrand, A0(n=0) and A2(n=0) terms 
f0_term0 = U*U*U * np.exp(-2.* U)\
        *2.*(1.+3.*a_1[0])*(1.+3.*a_2[0])
        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
f2_term0 =  U*U*U * np.exp(-2.* U)\
        *(1.-a_1[0])*(1.-a_2[0])
        #*((eiz_z[0]-eiz_w[0])/(eiz_z[0]-eiz_w[0]))\
Ft0_term0 = romb(f0_term0)
Ft2_term0 = romb(f2_term0)

# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
        # Integrand A0
        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
                *(2.*(1.+3.*a_1[j])*(1.+3.*a_2[j])*T*T*T*T\
                + 4.*(1.+2.*a_1[j]+2.*a_2[j]+3.*a_1[j]*a_2[j])*T*T \
                + 4.*(1.+a_1[j])*(1.+a_2[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        # Integrand A2                
        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
                /(T*T+1.)\
                *((T*T*T*T +4.*(T*T)+4.)*(1.-a_1[j])*(1.-a_2[j]))
                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
        Ft0 = romb(f0)
        Ft2 = romb(f2)
        #Integral = romb(D,dx=dqs,axis = 1)#intrgrate over dq
        
        A0[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
        #A0[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft0_term0
        A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta_1[j]*delta_2[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        #A2[i,0] = (1./2) * delta_1[0]*delta_2[0]*Ft2_term0
        A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0

        #A0[A0>1e6]= np.nan #nOTE: remove me later
        #A2[A2>1e6]= np.nan #nOTE: remove me later
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
    for k,theta in enumerate(thetas):
        G[i,k]=-(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*np.sin(thetas[k])*Ls[i]*Ls[i]*Ls[i]*Ls[i])*(sum_A0[i]+sum_A2[i]*np.cos(2.*thetas[k]))
        Tau[i,k]=-(np.pi*r_1*r_1*np.pi*r_2*r_2)/(2.*np.pi*Ls[i]*Ls[i]*Ls[i]*Ls[i])\
                *(sum_A0[i]+sum_A2[i]*(2.-np.cos(2.*thetas[k])))\
                *np.cos(thetas[k])/(np.sin(thetas[k])*np.sin(thetas[k]))
#np.savetxt('data/93alCG/eiz_x_CG.txt',eiz_x_1)
#np.savetxt('data/93alCG/eiz_z_CG.txt',eiz_z_1)
#np.savetxt('data/93alCG/eiz_x_93.txt',eiz_x_2)
#np.savetxt('data/93alCG/eiz_z_93.txt',eiz_z_2)
#np.savetxt('data/93alCG/eiz_i_AL.txt',eiz_w)
#np.savetxt('data/93al93/A0_93al93.txt',sum_A0)
#np.savetxt('data/93al93/A2_93al93.txt',sum_A2)
#np.savetxt('data/93al93/T_93al93.txt',Tau)
#np.savetxt('data/93al93/G_93al93.txt',G)

#np.savetxt('data/93alCG/A0_93alCG.txt',sum_A0)
#np.savetxt('data/93alCG/A2_93alCG.txt',sum_A2)
#np.savetxt('data/93alCG/T_93alCG.txt',Tau)
#np.savetxt('data/93alCG/G_93alCG.txt',G)


#pl.figure()
#pl.loglog(1e9*Ls,1e21*G[:,0])
#pl.loglog(1e9*Ls,1e21*G[:,1])
#pl.loglog(1e9*Ls,1e21*G[:,2])
#pl.title('Free energy vs separation')
#pl.show()
pl.figure()
pl.plot(thetas,1e21*G[0,:])
pl.plot(thetas,1e21*G[1,:])
pl.plot(thetas,1e21*G[2,:])
pl.title('Free energy vs mutual angle')
pl.show()
#pl.figure()
#pl.loglog(thetas,1e21*G[0,:])
#pl.loglog(thetas,1e21*G[1,:])
#pl.loglog(thetas,1e21*G[2,:])
#pl.title('Free energy vs mutual angle')
#pl.show()

pl.figure()
pl.plot(thetas,1e21*Tau[0,:])
pl.plot(thetas,1e21*Tau[1,:])
pl.plot(thetas,1e21*Tau[2,:])
pl.title('Torque vs mutual angle')
pl.show()
pl.figure()
pl.loglog(thetas,1e21*Tau[0,:])
pl.loglog(thetas,1e21*Tau[1,:])
pl.loglog(thetas,1e21*Tau[2,:])
pl.title('Torque vs mutual angle')
pl.show()

pl.figure()
pl.plot(1e9*Ls,1e21*Tau[:,0])
pl.plot(1e9*Ls,1e21*Tau[:,1])
pl.plot(1e9*Ls,1e21*Tau[:,2])
pl.title('Torque vs separation')
pl.show()
pl.figure()
pl.loglog(1e9*Ls,1e21*Tau[:,0])
pl.loglog(1e9*Ls,1e21*Tau[:,1])
pl.loglog(1e9*Ls,1e21*Tau[:,2])
pl.title('Torque vs separation')
pl.show()

pl.figure()
pl.plot(ns,eiz_z_2,'b-', label = r'$\epsilon^{93}_\hat{z}$')
pl.plot(ns,eiz_x_2,'b:', label = r'$\epsilon^{93}_\hat{x}$')
pl.plot(ns,eiz_z_1,'r-', label = r'$\epsilon^{CG10}_\hat{z}$')
pl.plot(ns,eiz_x_1,'r:', label = r'$\epsilon^{CG10}_\hat{x}$')
pl.plot(ns,eiz_w,'c-', label = r'$\epsilon^{Al2O3}$')
pl.title('LDS')
pl.legend(loc = 'best')
pl.axis([0,500,0,3.5])
#pl.savefig('plots/eiz_24_93_new_w.pdf')
pl.show()
pl.figure()
pl.loglog(ns,eiz_z_2,'b-', label = r'$\epsilon^{93}_\hat{z}$')
pl.loglog(ns,eiz_x_2,'b:', label = r'$\epsilon^{93}_\hat{x}$')
pl.loglog(ns,eiz_z_1,'r-', label = r'$\epsilon^{CG10}_\hat{z}$')
pl.loglog(ns,eiz_x_1,'r:', label = r'$\epsilon^{CG10}_\hat{x}$')
pl.loglog(ns,eiz_w,'c-', label = r'$\epsilon^{Al2O3}$')
pl.title('LDS')
pl.legend(loc = 'best')
#pl.savefig('plots/eiz_24_93_new_w.pdf')
pl.show()

pl.figure()
pl.plot(Ls,sum_A0/(sum_A0+sum_A2), 'b-', label = r'$A^{0}(\ell)/A$')
pl.plot(Ls,sum_A2/(sum_A0+sum_A2), 'r-', label = r'$A^{2}(\ell)/A$')
pl.plot(Ls,sum_A2-sum_A2, 'k:')
pl.plot(Ls,sum_A2/sum_A2, 'k:')
pl.legend(loc = 'best')
pl.title(r'Ratio of A0,A2 to A =(A0+A2)')
pl.ylabel(r'$\mathcal{A}^{(0)}/\mathcal{A}, \,\, \mathcal{A}^{(2)}/\mathcal{A}$')
pl.xlabel(r'$\ell\,\,\,[nm]$')
pl.show()

pl.figure()
pl.plot(1e9*Ls,1e21*sum_A0, 'b-', label = r'$A^{0}(\ell)$')
pl.plot(1e9*Ls,1e21*sum_A2, 'r-', label = r'$A^{2}(\ell)$')
pl.plot(1e9*Ls,(sum_A2-sum_A2), 'k:')
pl.legend(loc = 'best')
pl.title(r'Hamaker Coeffs A0,A2')
pl.ylabel(r'$\mathcal{A}^{(0)}, \,\, \mathcal{A}^{(2)}$')
pl.xlabel(r'$\ell\,\,\,[nm]$')
pl.show()

sys.exit()
pl.figure()
pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.title(r'93w93 Matsubara terms')
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/24w293_diff_delta2g2__A_vs_n.pdf')
pl.show()
#
pl.figure()
pl.plot(ns,(kbT/(32.)) * A2[0,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A2[1,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A2[2,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A2[3,:], 'b:', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/93_A_vs_n.pdf')
#pl.savefig('plots/24w93_nonlog_compare_delta2g2__A_terms_vs_n.pdf')
pl.show()
#
pl.figure()
pl.plot(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.plot(ns,(kbT/(32.)) * A2_org[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(kbT/(32.)) * A2_org[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.plot(ns,(kbT/(32.)) * A2_org[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.plot(ns,(kbT/(32.)) * A2_org[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/24w93_nonlog_A_terms_vs_n.pdf')
pl.show()

pl.figure()
pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/24w93_A_vs_n.pdf')
pl.show()

pl.figure()
pl.loglog(ns,(kbT/(32.)) * A2[ 0,:],'b-',label=r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[ 5,:],'g-',label=r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[5]))
pl.loglog(ns,(kbT/(32.)) * A2[10,:],'r-',label=r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[10]))
pl.loglog(ns,(kbT/(32.)) * A2[22,:],'c:')#, label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[24,:],'c-',label =r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[24]))
pl.loglog(ns,(kbT/(32.)) * A2[22,:],'c:')#, label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[47,:],'y:')#, label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[49,:],'y-',label =r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[49]))
pl.loglog(ns,(kbT/(32.)) * A2[51,:],'y:')#, label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(kbT/(32.)) * A2[99,:],'m-',label =r'$A^{(2)}(\ell=%2.1f\,nm)$'%(1e9*Ls[99]))
#
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(2)}_{n}$')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.legend(loc = 'best')
pl.axis([1.88,234.88,2.5e-36,1e-23])
#pl.savefig('plots/24w93_A_vs_n.pdf')
pl.show()
#
#
pl.figure()
pl.loglog(ns,A0[0,:]/sum_A0[0], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,A0[1,:]/sum_A0[1], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,A0[2,:]/sum_A0[2], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,A0[3,:]/sum_A0[3], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,A2[0,:]/sum_A2[0], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,A2[1,:]/sum_A2[1], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,A2[2,:]/sum_A2[2], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,A2[3,:]/sum_A2[3], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
pl.title(r'24w93 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
#pl.savefig('plots/24w93_relative_A_vs_n.pdf')
pl.show()

DA2_l = np.diff(A2, axis = 0)

DA2_n = np.diff(A2, axis = 1)

pl.figure()
pl.semilogx(ns,DA2_l[0,:],  'b-' ,linewidth = 1.0,label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[0]))
pl.semilogx(ns,DA2_l[9,:],  'b--',linewidth = 1.0,label =r'$\ell=%2.1f\,nm$' %(1e9*Ls[9]))
pl.semilogx(ns,DA2_l[24,:], 'b-.',linewidth = 1.0,label =r'$\ell=%2.1f\,nm$'%(1e9*Ls[24]))
pl.semilogx(ns,DA2_l[49,:], 'b:' ,linewidth = 1.0,label =r'$\ell=%2.1f\,nm$'%(1e9*Ls[49]))
pl.semilogx(ns,DA2_l[3,:]-DA2_l[3,:], 'k:')# , linewidth = 1.0, label = r'$\ell=%2.1f\,nm$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
pl.title(r'dA2/dl contribution to Matsubara sum for [5,1] and [26,0]')
pl.ylabel(r'$d\mathcal{A}^{(2)}_{n}/dl$')
pl.xlabel(r'$n$')
pl.show()

pl.figure()
pl.plot(ns,(A0[0,:]+A2[0,:]), 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.plot(ns,(A0[3,:]+A2[3,:]), 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.plot(ns,(A0[6,:]+A2[6,:]), 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[6]))
pl.plot(ns,(A0[9,:]+A2[9,:]), 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[9]))
pl.legend(loc = 'best')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.show()

pl.figure()
pl.loglog(ns,(A0[0,:]+A2[0,:])/(sum_A0[0]+sum_A2[0]), 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(A0[3,:]+A2[3,:])/(sum_A0[3]+sum_A2[3]), 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(A0[6,:]+A2[6,:])/(sum_A0[6]+sum_A2[6]), 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[6]))
pl.loglog(ns,(A0[9,:]+A2[9,:])/(sum_A0[9]+sum_A2[9]), 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[9]))
pl.legend(loc = 'best')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}, \,\, \mathcal{A}^{(2)}_{n}$')
pl.xlabel(r'$n$')
pl.show()
