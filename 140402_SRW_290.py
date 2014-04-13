#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
#eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
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

#filename1=sys.argv[1] # eiz_x
#filename2=sys.argv[2] # eiz_z
#usertag= 'SRW' #sys.argv[3]    
#
#mn = filename1.split('/')[-1].split('_')[2].replace('.txt','') 
## splits filename to get [m,n]
#
#titstr=usertag+' '+mn
#sys.exit()
#
eiz_x = np.loadtxt('data/eiz_x_290.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_290.txt') # LDS in parallel direction
#eiz_x = np.loadtxt(filename1) # LDS in perpendicular direction
#eiz_z = np.loadtxt(filename2) # LDS in parallel direction
##
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
##eiz_w[0] = 79.0

# Constants
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
#coeff = 0.159           # [eV]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]

# Matsubara frequencies
#ns = z/(0.159)#

ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(1e-9,100e-9,1e-9)  # separation distance between 2 cyclinders

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
pl.figure()
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
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
        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
        
        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0

        #A0[A0>1e6]= np.nan #NOTE: remove me later
        #A2[A2>1e6]= np.nan #NOTE: remove me later
    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
#pl.title(r'90w90 Matsubara terms')
#pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms')
pl.title(r'290w290 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
pl.xlabel(r'$N$')
#pl.savefig('plots/65_A_vs_n.pdf')
#pl.savefig('plots/90_A_vs_n.pdf')
#pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_A_vs_n.pdf')
pl.savefig('plots/290_A_vs_n.pdf')
pl.show()

print 'A0(separation) = ',sum_A0
print 'A2(separation) = ',sum_A2
print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]

#np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
#
#np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
#
#np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
#
#np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
#
np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)

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
	return 'plots/140322_%sw%s_HCs_%s_ret.pdf'%(cnt1,cnt2,orientation)

pl.figure()
pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[21],1e21*sum_A0[21]))
pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[21],1e21*sum_A2[21]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.legend(loc = 'best')

#pl.title(title('6','5','Log-log, perpendicular'))
#pl.title(title('9','0','perpendicular'))
#pl.title(title('9','1','perpendicular'))
#pl.title(title('9','3','perpendicular'))
pl.title(title('29','0','perpendicular'))
#
#pl.savefig(svfig('65','65','perpendicular'))
#pl.savefig(svfig('90','90','perpendicular'))
#pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','perpendicular'))
pl.savefig(svfig('290','290','perpendicular'))
pl.show()

pl.figure()
pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[21],1e21*sum_A0[21]))
pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[21],1e21*sum_A2[21]))
pl.xlabel(x_ax)
pl.ylabel(y_ax_per)
pl.title(title('6','5','perpendicular'))
pl.legend(loc = 'best')
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.legend(loc = 'best')
pl.title(title('6','5','Semi-log, perpendicular'))
#pl.title(title('9','0','perpendicular'))
#pl.title(title('9','1','perpendicular'))
#pl.title(title('9','3','perpendicular'))
#pl.title(title('29','0','perpendicular'))
#
pl.savefig(svfig('65','65','semilog_perpendicular'))
#pl.savefig(svfig('90','90','perpendicular'))
#pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','perpendicular'))
#pl.savefig(svfig('290','290','perpendicular'))
pl.show()

