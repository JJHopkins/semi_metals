#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

# Input dielectric response data
A0 = np.loadtxt('data/A0_p.txt') # LDS in perpendicular direction
A2 = np.loadtxt('data/A2_p.txt') # LDS in parallel direction
LS  = np.loadtxt('data/L.txt') # LDS in parallel direction
NS  = np.loadtxt('data/N.txt') # LDS in parallel direction
# Constants
c = 2.99e8              # [m/s]
coeff = 2.411e14        # [rad/s]
#coeff = 0.159           # [eV]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23 # [J]

ns = np.arange(0.,500.) 
zs = ns * coeff         

#Ls = np.arange(1e-9,101e-9,2e-9)  # separation distance between 2 cyclinders

from pylab import *
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##contourf(X,Y,h, 1000, cmap = hot())
#surf=ax.plot_surface(1e9*LS,NS,np.log(1e21*A0),rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
#CS = contour(1e9*LS,NS,np.log(1e21*A0), colors = 'k', linewidth = 0.5)
#pl.show()

X,Y = gradient(A0)

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
#pl.figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
cbar = pl.colorbar(surf)
cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
cbar.add_lines(CS)
#cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
ax.view_init(elev = 17, azim = 10)
savefig('plots/grad_A0_93_0.png')#, dpi = 300)
pl.show()

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
#pl.figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
cbar = pl.colorbar(surf)
cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
cbar.add_lines(CS)
#cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
ax.view_init(elev = 17, azim = 40)
savefig('plots/grad_A0_93_1.png')#, dpi = 300)
pl.show()

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
#pl.figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
cbar = pl.colorbar(surf)
cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
cbar.add_lines(CS)
#cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
ax.view_init(elev = 17, azim = 70)
savefig('plots/grad_A0_93_2.png')#, dpi = 300)
pl.show()

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
#pl.figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
cbar = pl.colorbar(surf)
cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
cbar.add_lines(CS)
#cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
ax.view_init(elev = 17, azim = 100)
savefig('plots/grad_A0_93_3.png')#, dpi = 300)
pl.show()

fig = pl.figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
#pl.figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
cbar = pl.colorbar(surf)
cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
cbar.add_lines(CS)
#cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
ax.view_init(elev = 17, azim = 130)
savefig('plots/grad_A0_93_4.png')#, dpi = 300)
pl.show()

#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
##pl.figure()
##contourf(X,Y,h, 1000, cmap = hot())
#surf=ax.plot_surface(1e9*LS,NS,Y,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
#CS = contour(1e9*LS,NS,Y, colors = 'k', linewidth = 0.5)
#pl.show()

pl.figure()
pl.contourf(1e9*LS,NS,X, cmap = cm.Blues)
#pl.quiver(1e9*LS[:,::10],NS[:,::10],X[:,::10], units = 'x',pivot='mid',width=0.02)
pl.quiver(1e9*LS,NS,X)#, units = 'x',pivot='mid')#,width=0.02)
#        np.log(A0[:,::10]), units = 'x',pivot='mid',width=0.02)
pl.show()
#    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
#    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
#pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
##pl.title(r'65w65 Matsubara terms')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
#pl.title(r'93w93 Matsubara terms with static term = 0')
##pl.title(r'290w290 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
#pl.xlabel(r'$N$')
##pl.savefig('plots/65_A_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
#pl.savefig('plots/93_n0_A_vs_n.eps')
##pl.savefig('plots/93_npk_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()
#
#print 'A0(separation) = ',sum_A0
#print 'A2(separation) = ',sum_A2
#print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
#print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
#
##np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
#np.savetxt('data/A0_93_n0_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_93_n0_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_93_n0_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_npk_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_npk_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_npk_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
#
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
#pl.title(title('9','3','A(n=0)=0, perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','n0_perpendicular'))
##pl.savefig(svfig('93','93','npk_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#pl.figure()
#pl.loglog(2.*Ls/c,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.xlabel(r'$2\ell/c$')
##pl.ylabel(r'(1/$\zeta_0$)2l/c')
#pl.legend(loc = 'best')
#
##pl.title(title('6','5','Log-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
#pl.title(title('9','3','A(n=0)=0, perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','deltaT_n0_perpendicular'))
##pl.savefig(svfig('93','93','npk_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#pl.figure()
#pl.plot(2.*np.sqrt(eiz_w[0])*Ls*zs[0]/c,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
#pl.ylabel(y_ax_per)
#pl.xlabel(r'$2\ell\sqrt{\varepsilon(i\zeta_{0})}/c$')
##pl.ylabel(r'(1/$\zeta_0$)2l/c')
#pl.legend(loc = 'best')
#
##pl.title(title('6','5','Log-log, perpendicular'))
##pl.title(title('9','0','perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
#pl.title(title('9','3','A(n=0)=0, perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','perpendicular'))
##pl.savefig(svfig('90','90','perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
#pl.savefig(svfig('93','93','p_n0_perpendicular'))
##pl.savefig(svfig('93','93','npk_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
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
#pl.title(title('9','3','A(n=0)=0,Semilog, perpendicular'))
##pl.title(title('9','1','perpendicular'))
##pl.title(title('9','3','perpendicular'))
##pl.title(title('29','0','perpendicular'))
##
##pl.savefig(svfig('65','65','semilog_perpendicular'))
##pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','npk_semilog_perpendicular'))
#pl.savefig(svfig('93','93','n0_semilog_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#
#
#
#

