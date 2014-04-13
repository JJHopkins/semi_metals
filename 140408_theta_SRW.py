#!/usr/bin/env python
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np
import matplotlib.pyplot as pl

lengths = np.loadtxt('data/Lengths_290_perpendicular_ret.txt')

A0_065 = np.loadtxt('data/A0_65_perpendicular_ret.txt')
A2_065 = np.loadtxt('data/A2_65_perpendicular_ret.txt')

A0_091 = np.loadtxt('data/A0_91_perpendicular_ret.txt')
A2_091 = np.loadtxt('data/A2_91_perpendicular_ret.txt')

A0_290 = np.loadtxt('data/A0_290_perpendicular_ret.txt')
A2_290 = np.loadtxt('data/A2_290_perpendicular_ret.txt')

thetas = np.linspace(0.0,np.pi/2,100)
X,Y = np.meshgrid(lengths,thetas)

A0_065_theta = A0_065 * np.cos(0.*Y) 
A0_091_theta = A0_091 * np.cos(0.*Y) 
A0_290_theta = A0_290 * np.cos(0.*Y) 

A2_065_theta = A2_065 * np.cos(2.*Y) 
A2_091_theta = A2_091 * np.cos(2.*Y) 
A2_290_theta = A2_290 * np.cos(2.*Y) 

##### A_2 PLOT S######
fig = figure()
#subplot(111, axisbg='darkslategray')
ax = fig.gca(projection = '3d')
#ax = fig.gca(projection = '3d', axisbg='darkslategray')

surf_065_0 = ax.plot_surface(1e9*X,Y,1e21*A0_065_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Blues, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()
surf_091_0 = ax.plot_surface(1e9*X,Y,1e21*A0_091_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Greens, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()
surf_290_0 = ax.plot_surface(1e9*X,Y,1e21*A0_290_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Reds, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()

CS_065_0 = contour(1e9*X,Y,1e21*A0_065_theta, colors = 'k', linewidth = 1.5)
CS_091_0 = contour(1e9*X,Y,1e21*A0_091_theta, colors = 'k', linewidth = 1.5)
CS_290_0 = contour(1e9*X,Y,1e21*A0_290_theta, colors = 'k', linewidth = 1.5)

cbar_065_0 = pl.colorbar(surf_065_0)
cbar_065_0.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[6,5]}\,\,\,[zJ]$', size = 14)
cbar_065_0.add_lines(CS_065_0)

cbar_091_0 = pl.colorbar(surf_091_0)
cbar_091_0.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,1]}\,\,\,[zJ]$', size = 14)
cbar_091_0.add_lines(CS_091_0)

cbar_290_0 = pl.colorbar(surf_290_0)
cbar_290_0.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[29,0]}\,\,\,[zJ]$', size = 14)
cbar_290_0.add_lines(CS_290_0)

cset_065_0 = ax.contour(1e9*X,Y,1e21*A0_065_theta,zdir='y',offset = -0.3,cmap = cm.Blues)
cset_091_0 = ax.contour(1e9*X,Y,1e21*A0_091_theta,zdir='y',offset = -0.3,cmap = cm.Greens)
cset_290_0 = ax.contour(1e9*X,Y,1e21*A0_290_theta,zdir='y',offset = -0.3,cmap = cm.Reds)
#man_loc = [(.1,.1),(.2,.2),(.3,.3),(.4,.4)]
yticks([0, pi/8, pi/6, pi/4, pi/3, pi/2],
        ['$0$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{\pi}{2}$'])
#clabel(CS, inline =1,fmt = '%1.1f', fontsize = 18,color = 'k', manual = man_loc)
ax.grid(on = True)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{angle}\,\,\,\theta\,\,[radians]$', size = 18)
ax.set_zlabel(r'$\mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{\mathcal{A}^{(0)}\, for \, [6,5],[9,1],\,and\,[29,0]\, in\,water}$',
        size = 21)
ax.view_init(elev = 10, azim = 65)
savefig('plots/A0_65_91_290.png')#, dpi = 300)
show()

##### A_2 PLOTS ######
fig = figure()
ax = fig.gca(projection = '3d')
surf_065 = ax.plot_surface(1e9*X,Y,1e21*A2_065_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Blues, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()
surf_091 = ax.plot_surface(1e9*X,Y,1e21*A2_091_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Greens, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()
surf_290 = ax.plot_surface(1e9*X,Y,1e21*A2_290_theta, rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Reds, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()

CS_065 = contour(1e9*X,Y,1e21*A2_065_theta, colors = 'k', linewidth = 1.5)
CS_091 = contour(1e9*X,Y,1e21*A2_091_theta, colors = 'k', linewidth = 1.5)
CS_290 = contour(1e9*X,Y,1e21*A2_290_theta, colors = 'k', linewidth = 1.5)

cbar_065 = pl.colorbar(surf_065)
cbar_065.ax.set_ylabel(r'$\mathcal{A}^{(2)}_{[6,5]}\,\,\,[zJ]$', size = 14)
cbar_065.add_lines(CS_065)

cbar_091 = pl.colorbar(surf_091)
cbar_091.ax.set_ylabel(r'$\mathcal{A}^{(2)}_{[9,1]}\,\,\,[zJ]$', size = 14)
cbar_091.add_lines(CS_091)

cbar_290 = pl.colorbar(surf_290)
cbar_290.ax.set_ylabel(r'$\mathcal{A}^{(2)}_{[29,0]}\,\,\,[zJ]$', size = 14)
cbar_290.add_lines(CS_290)

cset_065 = ax.contour(1e9*X,Y,1e21*A2_065_theta,zdir='x',offset = 110,cmap = cm.Blues)
cset_091 = ax.contour(1e9*X,Y,1e21*A2_091_theta,zdir='x',offset = 110,cmap = cm.Greens)
cset_290 = ax.contour(1e9*X,Y,1e21*A2_290_theta,zdir='x',offset = 110,cmap = cm.Reds)
#cset_065 = ax.contour(1e9*X,Y,1e21*A2_065_theta, zdir = 'y', offset =(1.1)*np.pi, cmap = cm.bone)# puts plot of max xi vs discrete r values at r=0 plane
#cset_065 = ax.contour(1e9*X,Y,1e21*A2_065_theta, zdir = 'z', offset = -1.1,       cmap = cm.bone)

#CS_065 = contour(1e9*X,Y,1e21*A2_065_theta, colors = 'k')
#man_loc = [(.1,.1),(.2,.2),(.3,.3),(.4,.4)]
#clabel(CS_065, inline =1,fmt = '%1.1f', fontsize = 18,color = 'k', manual = man_loc)
yticks([0, pi/8, pi/6, pi/4, pi/3, pi/2],
        ['$0$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{\pi}{2}$'])

ax.grid(on = True)

ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{angle}\,\,\,\theta\,\,[radians]$', size = 18)
ax.set_zlabel(r'$\mathcal{A}^{(2)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{\mathcal{A}^{(2)}\, for \, [6,5],[9,1],\,and\,[29,0]\, in\,water}$',
        size = 21)

#cbar_065 = pl.colorbar(surf_065)
#cbar_091 = pl.colorbar(surf_091)
#cbar_290 = pl.colorbar(surf_290)
#cbar_065.ax.set_ylabel(r'$\mathcal{A}^{(2)}_{[6,5]}\,\,\,[zJ]$', size = 14)
#cbar_065.add_lines(CS_065)
ax.view_init(elev = 12, azim = 154)
savefig('plots/A2_65_91_290.png')#, dpi = 300)
show()

##### A = (A_0 + A_2) PLOTS ######
fig = figure()
ax = fig.gca(projection = '3d')

surf065 = ax.plot_surface(1e9*X,Y,1e21*(A0_065_theta+A2_065_theta), rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Blues, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()

surf091 = ax.plot_surface(1e9*X,Y,1e21*(A0_091_theta+A2_091_theta), rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Greens, linewidth = 0.05, antialiased = True, shade = False)# True)#, cmap = hot()

surf290 = ax.plot_surface(1e9*X,Y,1e21*(A0_290_theta+A2_290_theta), rstride = 5, 
        cstride =5,alpha=0.7,cmap=cm.Reds, linewidth = 0.05, antialiased = True, shade = False)# True)True)#, cmap = hot()

cset065 = ax.contour(1e9*X,Y,1e21*(A0_065_theta+A2_065_theta),zdir='y',offset = -0.3,cmap = cm.Blues)
cset091 = ax.contour(1e9*X,Y,1e21*(A0_091_theta+A2_091_theta),zdir='y',offset = -0.3,cmap = cm.Greens)
cset290 = ax.contour(1e9*X,Y,1e21*(A0_290_theta+A2_290_theta),zdir='y',offset = -0.3,cmap = cm.Reds)
#cset = ax.contour(1e9*X,Y,1e21*A2_065_theta, zdir = 'y', offset =(1.1)*np.pi, cmap = cm.bone)# puts plot of max xi vs discrete r values at r=0 plane
#cset = ax.contour(1e9*X,Y,1e21*A2_065_theta, zdir = 'z', offset = -1.1,       cmap = cm.bone)

CS065 = contour(1e9*X,Y,1e21*(A0_065_theta+A2_065_theta), colors = 'k', linewidth = 1.5)
CS091 = contour(1e9*X,Y,1e21*(A0_091_theta+A2_091_theta), colors = 'k', linewidth = 1.5)
CS290 = contour(1e9*X,Y,1e21*(A0_290_theta+A2_290_theta), colors = 'k', linewidth = 1.5)

cbar065 = pl.colorbar(surf_065)
cbar065.ax.set_ylabel(r'$\mathcal{A}_{[6,5]}\,\,\,[zJ]$', size = 14)
cbar065.add_lines(CS_065)

cbar091 = pl.colorbar(surf_091)
cbar091.ax.set_ylabel(r'$\mathcal{A}_{[9,1]}\,\,\,[zJ]$', size = 14)
cbar091.add_lines(CS_091)

cbar290 = pl.colorbar(surf_290)
cbar290.ax.set_ylabel(r'$\mathcal{A}_{[29,0]}\,\,\,[zJ]$', size = 14)
cbar290.add_lines(CS_290)


#CS = contour(1e9*X,Y,1e21*A2_065_theta, colors = 'k')
#man_loc = [(.1,.1),(.2,.2),(.3,.3),(.4,.4)]
#clabel(CS, inline =1,fmt = '%1.1f', fontsize = 18,color = 'k', manual = man_loc)
yticks([0, pi/8, pi/6, pi/4, pi/3, pi/2],
        ['$0$', r'$\frac{\pi}{8}$', r'$\frac{\pi}{6}$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{3}$', r'$\frac{\pi}{2}$'])

ax.grid(on = True)
ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
ax.set_ylabel(r'$\rm{angle}\,\,\,\theta\,\,[radians]$', size = 18)
ax.set_zlabel(r'$(\mathcal{A}^{(0)}+\mathcal{A}^{(2)})\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
pl.title(r'$\rm{(\mathcal{A}^{(0)}+\mathcal{A}^{(2)})\, for \, [6,5],[9,1],\,and\,[29,0]\, in\,water}$',
        size = 21)

ax.view_init(elev = 10, azim = 65)
savefig('plots/A_tot_65_91_290.png')#, dpi = 300)
show()
