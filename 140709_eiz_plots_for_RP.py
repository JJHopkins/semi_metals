#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

eiz_x_24 = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z_24 = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
                                    
eiz_x_33 = np.loadtxt('data/eiz_x_33.txt') # LDS in perpendicular direction
eiz_z_33 = np.loadtxt('data/eiz_z_33.txt') # LDS in parallel direction
                                    
eiz_x_65 = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
eiz_z_65 = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction

eiz_x_90 = np.loadtxt('data/eiz_x_90.txt') # LDS in perpendicular direction
eiz_z_90 = np.loadtxt('data/eiz_z_90.txt') # LDS in parallel direction
        
eiz_x_91 = np.loadtxt('data/eiz_x_91.txt') # LDS in perpendicular direction
eiz_z_91 = np.loadtxt('data/eiz_z_91.txt') # LDS in parallel direction
         
eiz_x_93 = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
eiz_z_93 = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction

eiz_x_290= np.loadtxt('data/eiz_x_290.txt') # LDS in perpendicular direction
eiz_z_290= np.loadtxt('data/eiz_z_290.txt') # LDS in parallel direction

eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium

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

Ls = np.arange(1e-9,150e-9,1e-9)  # separation distance between 2 cyclinders
ls = np.arange(1e-9,1000e-9,1e-9)  # separation distance between 2 cyclinders
#pl.figure()
#pl.plot(zs,eiz_x_24, 'b-', label = r'$24 \epsilon_{x}$')
#pl.plot(zs,eiz_z_24, 'b--',label = r'$24 \epsilon_{z}$')
#pl.plot(zs,eiz_x_33, 'g-', label = r'$33 \epsilon_{x}$')
#pl.plot(zs,eiz_z_33, 'g--',label = r'$33 \epsilon_{z}$')
#pl.plot(zs,eiz_x_65, 'r-', label = r'$65 \epsilon_{x}$')
#pl.plot(zs,eiz_z_65, 'r--',label = r'$65 \epsilon_{z}$')
#pl.plot(zs,eiz_x_93, 'c-', label = r'$93\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_93, 'c--',label = r'$93\, \epsilon_{z},max= %2.2f$'%np.max(eiz_z_93))
#pl.plot(zs,eiz_w   , 'k-' ,label = r'$ \epsilon_{H2O}$')
#pl.legend(loc = 'best')
#pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
#pl.xlabel(r'$\omega_n \,\,\, [rad/sec]$', size = '24')
#pl.axis([0,1.22*1e17,0.9,10])
#pl.savefig('plots/eiz_24_33_65_93.pdf')
#pl.show()
#
#
#pl.figure()
#pl.loglog(zs,eiz_x_24, 'b-', label = r'$24$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_24, 'b:')#,label = r'$24\, \epsilon_{z}$')
#pl.loglog(zs,eiz_x_33, 'g-', label = r'$33$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_33, 'g:')#,label = r'$33\, \epsilon_{z}$')
#pl.loglog(zs,eiz_x_65, 'r-', label = r'$65$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_65, 'r:')#,label = r'$65\, \epsilon_{z}$')
#pl.loglog(zs,eiz_x_91, 'c-', label = r'$91$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_91, 'c:')#,label = r'$91\, \epsilon_{z}$')#,max= %2.2f$'%np.max(eiz_z_91))
#pl.loglog(zs,eiz_x_93, 'y-', label = r'$93$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_93, 'y:')#,label = r'$93\, \epsilon_{z},max= %2.2f$'%np.max(eiz_z_93))
#pl.loglog(zs,eiz_x_290,'m-', label = r'$290$')#\, \epsilon_{x}$')
#pl.loglog(zs,eiz_z_290,'m:')#,label = r'$290\, \epsilon_{z}$')#,max=%2.2f$'%np.max(eiz_z_290))
#pl.loglog(zs,eiz_w   , 'k-' ,label = r'$H_{2}O$')
#pl.legend(loc = 'best')
#pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
#pl.xlabel(r'$\omega_n \,\,\, [rad/sec]$', size = '24')
#pl.axis([0,5*1e17,0.9,12])
#pl.savefig('plots/loglog_eiz_24_33_65_93.pdf')
#pl.show()
#
#pl.figure()
#pl.plot(zs,eiz_x_24, 'b-', label = r'$24$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_24, 'b:')#,label = r'$24\, \epsilon_{z}$')
#pl.plot(zs,eiz_x_33, 'g-', label = r'$33$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_33, 'g:')#,label = r'$33\, \epsilon_{z}$')
#pl.plot(zs,eiz_x_65, 'r-', label = r'$65$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_65, 'r:')#,label = r'$65\, \epsilon_{z}$')
##plplotg(zs,eiz_x_90, 'c-', label = r'$90\, \epsilon_{x}$')
##plplotg(zs,eiz_z_90, 'c--')#,label = r'$90\, \epsilon_{z},max= %2.2f$'%np.max(eiz_z_90))
#pl.plot(zs,eiz_x_91, 'c-', label = r'$91$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_91, 'c:')#,label = r'$91\, \epsilon_{z}$')#,max= %2.2f$'%np.max(eiz_z_91))
#pl.plot(zs,eiz_x_93, 'y-', label = r'$93$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_93, 'y:')#,label = r'$93\, \epsilon_{z},max= %2.2f$'%np.max(eiz_z_93))
#pl.plot(zs,eiz_x_290,'m-', label = r'$290$')#\, \epsilon_{x}$')
#pl.plot(zs,eiz_z_290,'m:')#,label = r'$290\, \epsilon_{z}$')#,max=%2.2f$'%np.max(eiz_z_290))
#pl.plot(zs,eiz_w   , 'k-' ,label = r'$H_{2}O$')
#
#pl.legend(loc = 'best')
#pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
#pl.xlabel(r'$\omega_n \,\,\, [rad/sec]$', size = '24')
#pl.axis([0,1.22*1e17,0.9,10])
#pl.savefig('plots/all_eiz.pdf')
#pl.show()

A0_24 = np.loadtxt('data/A0_screen_24_sum.txt') # LDS in perpendicular direction
A2_24 = np.loadtxt('data/A2_screen_24_sum.txt') # LDS in parallel direction
A0_33 = np.loadtxt('data/A0_screen_33_sum.txt') # LDS in perpendicular direction
A2_33 = np.loadtxt('data/A2_screen_33_sum.txt') # LDS in parallel direction
A0_93 = np.loadtxt('data/A0_screen_93_sum.txt') # LDS in perpendicular direction
A2_93 = np.loadtxt('data/A2_screen_93_sum.txt') # LDS in parallel direction
A0_65 = np.loadtxt('data/65w65_srw_A0.txt') # LDS in perpendicular direction
A2_65 = np.loadtxt('data/65w65_srw_A2.txt') # LDS in parallel direction
A0_91 = np.loadtxt('data/91w91_srw_A0.txt') # LDS in perpendicular direction
A2_91 = np.loadtxt('data/91w91_srw_A2.txt') # LDS in parallel direction
A0_290= np.loadtxt('data/290w290_srw_A0.txt') # LDS in perpendicular direction
A2_290= np.loadtxt('data/290w290_srw_A2.txt') # LDS in parallel direction

pl.figure()
pl.loglog(1e9*Ls,A0_24 , 'b-', label = r'$24 $')#\, \epsilon_{x}$')
pl.loglog(1e9*Ls,A2_24 , 'b:')#,label = r'24 4\, \epsilon_{z}$')
pl.loglog(1e9*Ls,A0_33 , 'g-', label = r'$33 $')#\, \epsilon_{x}$')
pl.loglog(1e9*Ls,A2_33 , 'g:')#,label = r'33 3\, \epsilon_{z}$')
pl.loglog(1e9*Ls,A0_93 , 'r-', label = r'$93 $')#\, \epsilon_{x}$')
pl.loglog(1e9*Ls,A2_93 , 'r:')#,label = r'93 5\, \epsilon_{z}$')
pl.loglog(1e9*ls[:149],A0_65[:149] , 'c-', label = r'$65 $')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A2_65[:149] , 'c:')#,label = r'65 1\, \epsilon_{z}$')#,max= %2.2f$'%np.max(eiz_z_91))
pl.loglog(1e9*ls[:149],A0_91[:149] , 'y-', label = r'$91 $')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A2_91[:149] , 'y:')#,label = r'91 3\, \epsilon_{z},max= %2.2f$'%np.max(eiz_z_93))
pl.loglog(1e9*ls[:149],A0_290[:149],'m-', label = r'$290$')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A2_290[:149],'m:')#,label = r'$290290\, \epsilon_{z}$')#,max=%2.2f$'%np.max(eiz_z_290))
pl.legend(loc = 'best')
pl.ylabel(r'$A0, A2$', size = '24')
pl.xlabel(r'$\ell [nm]$', size = '24')
#pl.axis([0,5*1e17,0.9,12])
pl.savefig('plots/loglog_A0A2.pdf')
pl.show()

pl.figure()
pl.loglog(1e9*Ls,A0_24/A2_24 , 'b-', label = r'$24 $')#\, \epsilon_{x}$')
pl.loglog(1e9*Ls,A0_33/A2_33 , 'g-', label = r'$33 $')#\, \epsilon_{x}$')
pl.loglog(1e9*Ls,A0_93/A2_93 , 'r-', label = r'$93 $')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A0_65[:149]/A2_65[:149] , 'c-', label = r'$65 $')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A0_91[:149]/A2_91[:149] , 'y-', label = r'$91 $')#\, \epsilon_{x}$')
pl.loglog(1e9*ls[:149],A0_290[:149]/A2_290[:149],'m-', label = r'$290$')#\, \epsilon_{x}$')
pl.legend(loc = 'best')
pl.ylabel(r'$A0/A2$', size = '24')
pl.xlabel(r'$\ell [nm]$', size = '24')
#pl.axis([0,5*1e17,0.9,12])
pl.savefig('plots/loglog_ratio_A0A2.pdf')
pl.show()

pl.figure()
pl.plot(1e9*Ls,A0_24/A2_24 , 'b-', label = r'$24 $')#\, \epsilon_{x}$')
pl.plot(1e9*Ls,A0_33/A2_33 , 'g-', label = r'$33 $')#\, \epsilon_{x}$')
pl.plot(1e9*Ls,A0_93/A2_93 , 'r-', label = r'$93 $')#\, \epsilon_{x}$')
pl.plot(1e9*ls[:149],A0_65[:149]/A2_65[:149] , 'c-', label = r'$65 $')#\, \epsilon_{x}$')
pl.plot(1e9*ls[:149],A0_290[:149]/A2_290[:149],'m-', label = r'$290$')#\, \epsilon_{x}$')
pl.legend(loc = 'best')
pl.ylabel(r'$A0/A2$', size = '24')
pl.xlabel(r'$\ell [nm]$', size = '24')
#pl.axis([0,5*1e17,0.9,12])
pl.savefig('plots/ratio_A0A2.pdf')
pl.show()

pl.figure()
pl.plot(1e9*Ls,(A0_24-A2_24)/A0_24 , 'b-', label = r'$24 $')#\, \epsilon_{x}$')
pl.plot(1e9*Ls,(A0_33-A2_33)/A0_33 , 'g-', label = r'$33 $')#\, \epsilon_{x}$')
pl.plot(1e9*Ls,(A0_93-A2_93)/A0_93 , 'r-', label = r'$93 $')#\, \epsilon_{x}$')
pl.plot(1e9*ls[:149],( A0_65[:149]-A2_65[:149] )/ A0_65[:149], 'c-', label = r'$65 $')#\, \epsilon_{x}$')
pl.plot(1e9*ls[:149],(A0_290[:149]-A2_290[:149])/A0_290[:149],'m-', label = r'$290$')#\, \epsilon_{x}$')
pl.legend(loc = 'best')
pl.ylabel(r'$(A0-A2)/A0$', size = '24')
pl.xlabel(r'$\ell [nm]$', size = '24')
#pl.axis([0,5*1e17,0.9,12])
pl.savefig('plots/percent_A0A2.pdf')
pl.show()


eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
numbs = [1.0,2.0,3.0]
ns = np.arange(0.,500.) 
pl.figure()
for i,factor in enumerate(numbs):
    eiz_x[i] = numbs[i]*eiz_w
    eiz_z[i] = 2.*(eiz_x[i]-eiz_w)*eiz_w/(eiz_x[i]+eiz_w) + eiz_w
    pl.plot(ns,eiz_w)
    pl.plot(ns,eiz_x)
    pl.plot(ns,eiz_z)
pl.show()

pl.figure()
pl.plot(ns,eiz_z[0]/eiz_x[0])
pl.plot(ns,eiz_z[1]/eiz_x[1])
pl.plot(ns,eiz_z[2]/eiz_x[2])
pl.show()


eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_24.txt') # LDS of water, intervening medium
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.figure()
pl.plot(ns,eiz_z_24,'b-', label = '$\epsilon_{z}$')
pl.plot(ns,eiz_x,'g-',    label = '$\epsilon_{x}$')
pl.plot(ns,eiz_w,'k-',    label = '$\epsilon_{m}$')
pl.plot(ns,eiz_z, 'b:',   label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$')
pl.legend(loc = 'best')
pl.show()

pl.figure()
pl.loglog(ns,eiz_z_24,'b-' ,label = '\epsilon_{z}')                         
pl.loglog(ns,eiz_x,'g-'    ,label = '\epsilon_{x}')                       
pl.loglog(ns,eiz_w,'k-'    ,label = '\epsilon_{m}')                      
pl.loglog(ns,eiz_z, 'b:'   ,label = '\epsilon_{z}(\mathcal{A}^{(2)}=0)') 
pl.legend(loc = 'best')
pl.show()


eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_93.txt') # LDS of water, intervening medium
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.figure()
pl.plot(ns,eiz_z_24,'b-', label = '$\epsilon_{z}$')
pl.plot(ns,eiz_x,'g-',    label = '$\epsilon_{x}$')
pl.plot(ns,eiz_w,'k-',    label = '$\epsilon_{m}$')
pl.plot(ns,eiz_z, 'b:',   label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$')
pl.legend(loc = 'best')
pl.show()

pl.figure()
pl.loglog(ns,eiz_z_24,'b-' ,label = '\epsilon_{z}')                         
pl.loglog(ns,eiz_x,'g-'    ,label = '\epsilon_{x}')                       
pl.loglog(ns,eiz_w,'k-'    ,label = '\epsilon_{m}')                      
pl.loglog(ns,eiz_z, 'b:'   ,label = '\epsilon_{z}(\mathcal{A}^{(2)}=0)') 
pl.legend(loc = 'best')
pl.show()


eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_65.txt') # LDS of water, intervening medium
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.figure()
pl.plot(ns,eiz_z_24,'b-', label = '$\epsilon_{z}$')
pl.plot(ns,eiz_x,'g-',    label = '$\epsilon_{x}$')
pl.plot(ns,eiz_w,'k-',    label = '$\epsilon_{m}$')
pl.plot(ns,eiz_z, 'b:',   label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$')
pl.legend(loc = 'best')
pl.show()

pl.figure()
pl.loglog(ns,eiz_z_24,'b-' ,label = '\epsilon_{z}')                         
pl.loglog(ns,eiz_x,'g-'    ,label = '\epsilon_{x}')                       
pl.loglog(ns,eiz_w,'k-'    ,label = '\epsilon_{m}')                      
pl.loglog(ns,eiz_z, 'b:'   ,label = '\epsilon_{z}(\mathcal{A}^{(2)}=0)') 
pl.legend(loc = 'best')
pl.show()



pl.figure()

#c = loadtxt('output/130807_3D_A_dep.txt', unpack=True, usecols = [0])
eiz_x = np.loadtxt('data/eiz_x_24.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_24.txt') # LDS of water, intervening medium
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.loglog(ns,eiz_z_24,'b-' ,label = '$24\,\epsilon_{z}(\mathcal{A}^{(2)}=0.0045)$')                         
pl.loglog(ns,eiz_z  , 'b:'   ,label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$') 

eiz_x = np.loadtxt('data/eiz_x_33.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_33.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.loglog(ns,eiz_z_24,'g-' ,label = '$33\,\epsilon_{z}(\mathcal{A}^{(2)}=0.84)$')                         
pl.loglog(ns,eiz_z  , 'g:'   ,label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$') 

eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS of water, intervening medium
eiz_z_24 = np.loadtxt('data/eiz_z_93.txt') # LDS of water, intervening medium
eiz_z = 2.*(eiz_x-eiz_w)*eiz_w/(eiz_x+eiz_w) + eiz_w

pl.loglog(ns,eiz_z_24,'r-' ,label = '$93\,\epsilon_{z}(\mathcal{A}^{(2)}=3.23)$')                         
pl.loglog(ns,eiz_z  , 'r:'   ,label = '$\epsilon_{z}(\mathcal{A}^{(2)}=0)$') 

x_93,eiz_DD_93 = np.loadtxt('data/SWCNT_93_z.txt', unpack=True, usecols = [0,1])
x_24, eiz_DD_24 = np.loadtxt('data/SWCNT_24_z.txt', unpack=True, usecols =[0,1])

pl.loglog(x_24/0.159,eiz_DD_24,'b--' ,label = '$24\,\epsilon_{z}(Gecko)$')                         
pl.loglog(x_93/0.159,eiz_DD_93,'r--'   ,label = '$93\,\epsilon_{z}(Gecko)$') 

pl.legend(loc = 'best')
#pl.axis([0,502,0,20])
pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
pl.xlabel(r'$n$', size = '24')
pl.savefig(plots/'eiz_CNT.pdf')
pl.show()



eV_DD, eiz_wDD  = np.loadtxt('data/water_fine_DD_140708.txt', unpack=True, usecols = [0,1])
x    ,        y = np.loadtxt('data/water_finer_DD_140708.txt', unpack=True, usecols =[0,1])

pl.loglog(x_24/0.159,eiz_DD_24,'b--' ,label = '$24\,\epsilon_{z}(Gecko)$')                         
pl.loglog(x_93/0.159,eiz_DD_93,'r--'   ,label = '$93\,\epsilon_{z}(Gecko)$') 

pl.legend(loc = 'best')
#pl.axis([0,502,0,20])
pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
pl.xlabel(r'$n$', size = '24')
pl.savefig(plots/'eiz_CNT.pdf')
pl.show()
