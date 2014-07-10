#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
#Jx,Jy = np.loadtxt('data/im_Jcv.csv', delimiter = ',',unpack=True, usecols = [0,1]) # LDS in perpendicular direction
Jx,Jy = np.loadtxt('data/sio2-amorp-t.csv', delimiter = ',',unpack=True, usecols = [0,1]) # LDS in perpendicular direction
#ex,ey = np.loadtxt('data/a-sio2-eps2.txt',unpack=True, usecols = [0,1])
ex,ey = np.loadtxt('data/sio2-amorp-tKL.csv', delimiter = ',',unpack=True, usecols = [0,1]) # LDS in perpendicular direction
#eiz_z = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction
#
#inJ = 1./np.sqrt(Jy[296:])
inJ = 1./Jy
#inJ[:296] = 0
inJ[:370] = 0
multJ = 8.289*1e6*Jy

#eps2 = 8.289*1e6*Jy/(Jx*Jx)
eps2 = Jy/(Jx*Jx)

pl.figure()
pl.plot(Jx,eps2)
pl.plot(Jx,1e-2*Jy)
#pl.plot(ex,ey)
pl.title('eps2')
pl.show()

#pl.figure()
#pl.plot(Jx,Jy)
#pl.plot(Jx,8.289*1e-6*Jy)
#pl.plot(ex,ey)
#pl.title('Jcv')
#pl.show()
#
#pl.figure()
#pl.plot(Jx,eps2)
#pl.plot(ex,ey)
#pl.title('eps2')
#pl.show()
#
#pl.figure()
#pl.plot(Jx,np.sqrt(Jy))
#pl.plot(Jx,np.sqrt(8.289*1e6*Jy))
#pl.plot(Jx,8.289*1e-6*np.sqrt(Jy))
#pl.plot(ex,ey)
#pl.title('sqrt Jcv')
#pl.show()
#
#pl.figure()
#pl.plot(Jx,inJ)
#pl.plot(ex,ey)
#pl.title('1/Jcv')
##pl.axis([0,50,0,55])
#pl.show()
#
#pl.figure()
#pl.plot(Jx,1./multJ)
#pl.plot(ex,ey)
#pl.title('1e6*Jcv')
#pl.axis([0,50,0,25])
#pl.show()

eiz_x_24 = np.loadtxt('data/eiz_x_24.txt') # LDS in perpendicular direction
eiz_z_24 = np.loadtxt('data/eiz_z_24.txt') # LDS in parallel direction

eiz_x_93 = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
eiz_z_93 = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction

eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium


pl.figure()
pl.plot(zs,eiz_x_24, 'g-', label = r'$24 \epsilon_{x}$')
pl.plot(zs,eiz_z_24, 'g--', label = r'$24 \epsilon_{z}$')

pl.plot(zs,eiz_x_93, 'r-', label =  r'$93\, \epsilon_{x}$')
pl.plot(zs,eiz_z_93, 'r--', label = r'$93\, \epsilon_{z}$')
pl.plot(zs,eiz_z_93, ' ', label = r'$max \epsilon_{z} = %2.2f$'%np.max(eiz_z_93))

pl.plot(zs,eiz_w, 'b-', label = r'$93 \epsilon_{x}$')
pl.legend(loc = 'best')
pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
pl.xlabel(r'$\omega_n \,\,\, [rad/sec]$', size = '24')
pl.savefig('compare_93_24.pdf')
pl.show()


pl.figure()
pl.loglog(zs,eiz_x_24, 'g-',  label = r'$24\, \epsilon_{x}$')
pl.loglog(zs,eiz_z_24, 'g--', label = r'$24\, \epsilon_{z}$')
pl.loglog(zs,eiz_x_93, 'r-',  label =  r'$93\, \epsilon_{x}$')
pl.loglog(zs,eiz_z_93, 'r--', label = r'$93\, \epsilon_{z}$')
#pl.loglog(zs,eiz_z_93, ' ',   label = r'$max[\epsilon_{z}] = %2.2f$'%np.max(eiz_z_93))
pl.loglog(zs,eiz_w,     'b-',  label = r'$93 \epsilon_{x}$')
pl.legend(loc = 'best')
pl.ylabel(r'$\epsilon(i\omega)$', size = '24')
pl.xlabel(r'$\omega_n \,\,\, [rad/sec]$', size = '24')
pl.savefig('loglog_compare_93_24.pdf')
pl.show()
