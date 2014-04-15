#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

#
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium
eiz_x = np.loadtxt('data/eiz_x_93.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_93.txt') # LDS in parallel direction
Ls = np.arange(1e-9,500e-9,50e-9)  # separation distance between 2 cyclinders

A0 = np.loadtxt('data/A0_93_perpendicular_ret.txt')
A2 = np.loadtxt('data/A2_93_perpendicular_ret.txt')

A0_n0 = np.loadtxt('data/A0_93_n0_perpendicular_ret.txt')
A2_n0 = np.loadtxt('data/A2_93_n0_perpendicular_ret.txt')
#

pl.figure()
pl.loglog(Ls,A0,'b-')
pl.loglog(Ls,A0_n0,'g-')
pl.show()
pl.figure()
pl.loglog(Ls,A2,'b--')
pl.loglog(Ls,A2_n0,'g--')
pl.show()


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
#	return 'plots/140322_%sw%s_HCs_%s_ret.eps'%(cnt1,cnt2,orientation)
#
#pl.figure()
#pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
#pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.legend(loc = 'best')
#
##pl.title(title('6','5','Log-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
#pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','perpendicular'))
##pl.savefig(svfig('93','93','npk_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#pl.figure()
#pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'      %(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
#pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
#pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.legend(loc = 'best')
##pl.title(title('6','5','Semi-log, perpendicular'))
##pl.title(title('9','3','(no first peak eps2_z) Semi-log, perpendicular'))
#pl.title(title('9','3','Semi-log, perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','semilog_perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','npk_semilog_perpendicular'))
#pl.savefig(svfig('93','93','semilog_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#
#
