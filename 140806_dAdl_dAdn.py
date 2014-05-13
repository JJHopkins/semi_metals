#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import *
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
#coeff = 0.159           # [eV]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]

# Matsubara frequencies
#ns = z/(0.159)#


x_labels =r'$\rm{separation}\,\,\,\ell\,\,[nm]$'#,   size = 18      #) r'$\,\ell\,\,\,\rm{[nm]}$'
y_labels = r'$\rm{n^{th}}\,Matsubra\,term$'     #, size = 18         #      r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
z_labels = r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$'#,size = 18 ,rotation ='horizontal'

def title(cnt_m,cnt_n):
	return r'$\rm{[%s,%s]\,Gradient\,\mathcal{A}^{(0)}\,terms:\,perpendicular\,in\,water,\,retarded}$'%(cnt_m,cnt_n)

def svfig(mn):
	return 'plots/grad_A0_%s.png'%(mn)

def proj_val(gradmin):
    return min(gradmin) - (0.1)*gradmin

def plsvfig(LLS,NNS,X,projval,x_label,y_label,z_label,Title_fnc,save_fnc):
    fig = pl.figure()
    ax = fig.gca(projection = '3d')
    surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
    CS = contour(1e9*LLS,NNS,X, colors = 'k', linewidth = 0.5)
    cbar = pl.colorbar(surf)
    #cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
    cbar.add_lines(CS)
    cset = ax.contour(1e9*LLS,NNS,X,zdir='z',offset = projval ,cmap = cm.Blues)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    pl.title(Title_fnc)
    #pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
    ax.view_init(elev = 17, azim = 150)
    savefig(save_fnc)#'plots/grad_A0_93_0.png')#, dpi = 300)
    #pl.show()
    return 0

def cplsvfig(LLS,NNS,X,projval,x_label,y_label,z_label,Title_fnc,save_fnc):
    surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.bone,linewidth = 0.05, antialiased = True, shade = False)
    CS = contour(1e9*LLS,NNS,X, colors = 'k', linewidth = 0.5)
    cbar = pl.colorbar(surf)
    #cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
    cbar.add_lines(CS)
    cset = ax.contour(1e9*LLS,NNS,X,zdir='z',offset = projval ,cmap = cm.Blues)
    return 0

ns = np.arange(0.,500.) 
nns =np.arange(0,499.) 
Ls = np.arange(1e-9,450e-9,1e-9)  # separation distance between 2 cyclinders
LLs= Ls[:448]
#ns = np.arange(0.,300.) 
#Ls = np.arange(1e-9,500e-9,10e-9)  # separation distance between 2 cyclinders
NS,LS = np.meshgrid(ns,Ls)
NNs,LLLs = np.meshgrid(nns,Ls)

A0n65  = np.loadtxt('data/A0_n_65.txt') # LDS in perpendicular direction
A0sum65= np.loadtxt('data/A0_65_sum.txt') # LDS in perpendicular direction
A0n93  = np.loadtxt('data/A0_n_93.txt') # LDS in perpendicular direction
A0sum93= np.loadtxt('data/A0_93_sum.txt') # LDS in perpendicular direction

dAdl65 = (kbT/32)*np.diff(A0n65,axis = 0)
dAdn65 = (kbT/32)*np.diff(A0n65,axis = 1)

dAdl93 = (kbT/32)*np.diff(A0n93,axis = 0)
dAdn93 = (kbT/32)*np.diff(A0n93,axis = 1)

dsumAdl65 = np.diff(A0sum65,axis = 0)
dsumAdl93 = np.diff(A0sum93,axis = 0)

#A0n65_trunc050 = A0n65[:,:100]
#A0sum65_trunc050 = np.sum(A0n65_trunc050, axis = 1)

A0n65_trunc150 = A0n65[:,:150]
A0sum65_trunc150 = np.sum(A0n65_trunc150, axis = 1)

A0n65_trunc300 = A0n65[:,:300]
A0sum65_trunc300 = np.sum(A0n65_trunc300, axis = 1) 

#dAsum65_trunc050dl = (kbT/32)*np.diff(A0sum65_trunc050)
dAsum65_trunc150dl = (kbT/32)*np.diff(A0sum65_trunc150)
dAsum65_trunc300dl = (kbT/32)*np.diff(A0sum65_trunc300)

pl.figure()
pl.plot(1e9*LLs,1e21*dsumAdl65, 'b-', label = r'$n_{max} = 500$')
pl.plot(1e9*LLs,1e21*dAsum65_trunc300dl, 'g-', label = r'$n_{max} = 300$')
pl.plot(1e9*LLs,1e21*dAsum65_trunc150dl, 'r-', label = r'$n_{max} = 150$')
#pl.plot(LLs,dAsum65_trunc050dl, 'c-', label = r'$n_{max} = 100$')
pl.legend(loc = 'best')      
#pl.title(r'65w65 p(n,l)')
pl.ylabel(r'd$\cal{A}^{(0)}$/d$\ell$   [zJ]')
pl.xlabel(r'$n$   [nm]')
pl.axis([0,1e9*1.5e-7,-2.0,0.0])
pl.savefig('plots/65_dsumAdl_vs_n.png')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
pl.show()                    

pl.figure()
pl.plot(nns,1e21*dAdn65[ 0,:]  , 'b-', label = r'$\ell$= 1 to 11 nm')
pl.plot(nns,1e21*dAdn65[ 5,:]  , 'b-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))
pl.plot(nns,1e21*dAdn65[10,:]  , 'b-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))

pl.plot(nns,1e21*dAdn65[50,:]  , 'g-', label = r'$\ell$= 50 to 61 nm')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
pl.plot(nns,1e21*dAdn65[55,:]  , 'g-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
pl.plot(nns,1e21*dAdn65[60,:]  , 'g-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))

pl.plot(nns,1e21*dAdn65[100,:]  , 'r-', label = r'$\ell$= 100 to 111 nm')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
pl.plot(nns,1e21*dAdn65[105,:]  , 'r-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
pl.plot(nns,1e21*dAdn65[110,:]  , 'r-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))

pl.plot(nns,1e21*dAdn65[150,:]  , 'c-', label = r'$\ell$= 150 to 161 nm')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
pl.plot(nns,1e21*dAdn65[155,:]  , 'c-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
pl.plot(nns,1e21*dAdn65[160,:]  , 'c-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))

pl.plot(nns,1e21*dAdn65[200,:] , 'y-', label = r'$\ell$= 200 to 211 nm')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[50]))
pl.plot(nns,1e21*dAdn65[205,:] , 'y-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[50]))
pl.plot(nns,1e21*dAdn65[210,:] , 'y-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[50]))

pl.plot(nns,1e21*dAdn65[250,:] , 'm-', label = r'$\ell$= 250 to 261 nm')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[100]))
pl.plot(nns,1e21*dAdn65[255,:] , 'm-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[100]))
pl.plot(nns,1e21*dAdn65[260,:] , 'm-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[100]))
#pl.plot(nns,dAdn65[260,:] , 'k-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[200]))
#pl.plot(nns,dAdn65[255,:] , 'k-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[200]))
#pl.plot(nns,dAdn65[260,:] , 'k-')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[200]))
#pl.plot(NNs[200:205,:],dAdn65[200:205,:], 'g-.')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[300]))
#pl.plot(NNs[300:305,:],dAdn65[300:305,:], 'r-.')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[400]))
#pl.plot(nns,p[25,:], 'c-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[25]))
#pl.plot(nns,p[30,:], 'b-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
#pl.plot(nns,p[35,:], 'g-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[35]))
#pl.plot(nns,p[40,:], 'r-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
#pl.plot(nns,p[45,:], 'c-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[45]))
#pl.plot(nns, dAdn65[10,:]/dAdn65[10,:], 'k-')# , label = r'$p(l=%2.1f)$'%(Ls[45]))
pl.axis([0,37,-0.8,0.3])
pl.legend(loc = 'best')      
#pl.title(r'65w65 p(n,l)')
pl.ylabel(r'd$\cal{A}_{n}^{(0)}$/dn    [zJ]')
pl.xlabel(r'$n$')
pl.savefig('plots/65_dAdn_vs_n.png')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
pl.show()                    

#pl.figure()
##pl.plot(NNs[0:5,:]    ,dAdn65[0:5,:]    , 'b-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))
##pl.plot(NNs[20:25,:]  ,dAdn65[20:25,:]  , 'g-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
##pl.plot(NNs[50:55,:]  ,dAdn65[50:55,:]  , 'r-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
#pl.plot(NNs[90:95,:]  ,dAdn65[90:95,:]  , 'c-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
#pl.plot(NNs[95:100,:] ,dAdn65[95:100,:] , 'y-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[50]))
#pl.plot(NNs[100:105,:],dAdn65[100:105,:], 'm-' )#, label = r'$p(l=%2.1f)$'%(1e9*Ls[100]))
#pl.plot(NNs[105:110,:],dAdn65[105:110,:], 'b-.')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[200]))
##pl.plot(NNs[200:205,:],dAdn65[200:205,:], 'g-.')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[300]))
##pl.plot(NNs[300:305,:],dAdn65[300:305,:], 'r-.')#, label = r'$p(l=%2.1f)$'%(1e9*Ls[400]))
##pl.plot(nns,p[25,:], 'c-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[25]))
##pl.plot(nns,p[30,:], 'b-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
##pl.plot(nns,p[35,:], 'g-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[35]))
##pl.plot(nns,p[40,:], 'r-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
##pl.plot(nns,p[45,:], 'c-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[45]))
#pl.plot(nns, dAdn65[10,:]/dAdn65[10,:], 'k-')# , label = r'$p(l=%2.1f)$'%(Ls[45]))
#pl.legend(loc = 'best')      
##pl.title(r'65w65 p(n,l)')
#pl.ylabel(r'd$\cal{A}_{n}^{(0)}$/dn')
#pl.xlabel(r'$n$')
#pl.savefig('plots/65_dAdn_vs_n.png')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()                    



## Input dielectric response data
#A0_65  = np.loadtxt('data/A0_65_p.txt') # LDS in perpendicular direction
#A0_91  = np.loadtxt('data/A0_91_p.txt') # LDS in perpendicular direction
#A0_93  = np.loadtxt('data/A0_93_p.txt') # LDS in perpendicular direction
#A0_290 = np.loadtxt('data/A0_290_p.txt') # LDS in perpendicular direction
#
#X_65,Y_65   = gradient(A0_65)
#X_91,Y_91   = gradient(A0_91)
#X_93,Y_93   = gradient(A0_93)
#X_290,Y_290 = gradient(A0_290)
#
#plsvfig(LS,NS,X_65, proj_val(X_65[:, 25] ), x_labels,y_labels,z_labels,title('6','5'), svfig('65'))
#plsvfig(LS,NS,X_91, proj_val(X_91[:, 25] ), x_labels,y_labels,z_labels,title('9','1'), svfig('91'))
#plsvfig(LS,NS,X_93, proj_val(X_93[:, 25] ), x_labels,y_labels,z_labels,title('9','3'), svfig('93'))
#plsvfig(LS,NS,X_290,proj_val(X_290[:,25]),x_labels,y_labels,z_labels,title('29','0'),svfig('290'))
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#cplsvfig(LS,NS,X_65, proj_val(X_65[:, 25] ), x_labels,y_labels,z_labels,title('6','5'), svfig('65'))
#cplsvfig(LS,NS,X_91, proj_val(X_91[:, 25] ), x_labels,y_labels,z_labels,title('9','1'), svfig('91'))
#cplsvfig(LS,NS,X_93, proj_val(X_93[:, 25] ), x_labels,y_labels,z_labels,title('9','3'), svfig('93'))
#cplsvfig(LS,NS,X_290,proj_val(X_290[:,25]),x_labels,y_labels,z_labels,title('29','0'),svfig('290'))
##ax.set_xlabel(x_label)
##ax.set_ylabel(y_label)
##ax.set_zlabel(z_label)
##pl.title(Title_fnc)
###pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
#ax.view_init(elev = 17, azim = 150)
#savefig('plots/combo_contour.png')#'plots/grad_A0_93_0.png')#, dpi = 300)
##pl.show()
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#surf=ax.plot_surface(1e9*LS,NS,X_65,rstride=1,cstride=1, alpha=0.5,cmap=cm.Blues,linewidth = 0.05, antialiased = True, shade = False)
#surf=ax.plot_surface(1e9*LS,NS,X_91,rstride=1,cstride=1, alpha=0.5,cmap=cm.Greens,linewidth = 0.05, antialiased = True, shade = False)
#surf=ax.plot_surface(1e9*LS,NS,X_93,rstride=1,cstride=1, alpha=0.5,cmap=cm.Reds,linewidth = 0.05, antialiased = True, shade = False)
#surf=ax.plot_surface(1e9*LS,NS,X_290,rstride=1,cstride=1,alpha=0.5,cmap=cm.hot,linewidth = 0.05, antialiased = True, shade = False)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
#ax.view_init(elev = 17, azim = 150)
##fig = pl.figure()
##ax = fig.gca(projection = '3d')
##surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
##ax.view_init(elev = 17, azim = 10)
#savefig('plots/grad_A0_combo.png')#, dpi = 300)
##show()
##
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#contour(1e9*LS,NS,X_65, rstride=1,cstride=1,alpha=0.5,cmap=cm.Blues,  linewidth = 0.05, antialiased = True, shade = False)
#contour(1e9*LS,NS,X_91, rstride=1,cstride=1,alpha=0.5,cmap=cm.Greens, linewidth = 0.05, antialiased = True, shade = False)
#contour(1e9*LS,NS,X_93, rstride=1,cstride=1,alpha=0.5,cmap=cm.Reds,   linewidth = 0.05, antialiased = True, shade = False)
#contour(1e9*LS,NS,X_290,rstride=1,cstride=1,alpha=0.5,cmap=cm.Purples,linewidth = 0.05, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}}$',size = 21)
##ax.view_init(elev = 17, azim = 40)
#savefig('plots/A0_contours.png')#, dpi = 300)
##pl.show()
##
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#ax.plot_wireframe(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="blue")
##ax.plot_wireframe(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="green")
#ax.plot_wireframe(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.4,linewidth=0.4,color="red")
##ax.plot_wireframe(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="magenta")
#contour(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.5, colors = 'b')#,alpha=0.7,color="blue"    ,linewidth = 0.6, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.7,color="green"   ,linewidth = 0.6, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.9, colors = 'r')#,alpha=0.7,color="red"     ,linewidth = 0.6, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.7,color="magenta" ,linewidth = 0.6, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$-\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{-Gradient\, \mathcal{A}^{(0)}}$',size = 21)
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
###ax.view_init(elev = 16, azim = 175)
#savefig('plots/A0_wire2_0.png')#, dpi = 300)
#pl.show()
##
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.plot_wireframe(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="blue")
#ax.plot_wireframe(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="green")
##ax.plot_wireframe(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.4,linewidth=0.4,color="red")
#ax.plot_wireframe(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="magenta")
##contour(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.7,color="blue"    ,linewidth = 0.7, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.7, colors = 'g')#,color="green"   ,linewidth = 0.7, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.7,color="red"     ,linewidth = 0.7, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.9, colors = 'm')#,color="magenta" ,linewidth = 0.7, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$-\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{-Gradient\, \mathcal{A}^{(0)}}$',size = 21)
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
###ax.view_init(elev = 16, azim = 200)
#savefig('plots/A0_wire2_1.png')#, dpi = 300)
#pl.show()
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
#ax.plot_wireframe(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="blue")
##ax.plot_wireframe(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="green")
##ax.plot_wireframe(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.4,linewidth=0.4,color="red")
#ax.plot_wireframe(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="magenta")
#contour(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.5, colors = 'b')#,color='b' ,linewidth = 0.5, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.7,color='g' ,linewidth = 0.5, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.7,color='r' ,linewidth = 0.5, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.9, colors = 'm')#,color='m' ,linewidth = 0.5, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$-\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{-Gradient\, \mathcal{A}^{(0)}}$',size = 21)
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
###ax.view_init(elev = 16, azim = 225)
#savefig('plots/A0_wire2_2.png')#, dpi = 300)
#pl.show()
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.plot_wireframe(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="blue")
#ax.plot_wireframe(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="green")
#ax.plot_wireframe(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.4,linewidth=0.4,color="red")
##ax.plot_wireframe(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.3,linewidth=0.3,color="magenta")
##contour(1e9*LS,NS,-X_65, rstride=1,cstride=1,alpha=0.7,color='b' ,linewidth = 0.5, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_91, rstride=1,cstride=1,alpha=0.7, colors = 'g')#,color='g' ,linewidth = 0.5, antialiased = True, shade = False)
#contour(1e9*LS,NS,-X_93, rstride=1,cstride=1,alpha=0.9, colors = 'r')#,color='r' ,linewidth = 0.5, antialiased = True, shade = False)
##contour(1e9*LS,NS,-X_290,rstride=1,cstride=1,alpha=0.7,color='m' ,linewidth = 0.5, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
##cbar = pl.colorbar(surf)
##cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
##cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
##cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$-\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{-Gradient\, \mathcal{A}^{(0)}}$',size = 21)
##pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
###ax.view_init(elev = 16, azim = 225)
#savefig('plots/A0_wire2_3.png')#, dpi = 300)
#pl.show()
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
##pl.figure()
##contourf(X,Y,h, 1000, cmap = hot())
#surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
#CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
#cbar = pl.colorbar(surf)
#cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
#cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
#cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
#ax.view_init(elev = 17, azim = 70)
#savefig('plots/grad_A0_93_2.png')#, dpi = 300)
#pl.show()
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
##pl.figure()
##contourf(X,Y,h, 1000, cmap = hot())
#surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
#CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
#cbar = pl.colorbar(surf)
#cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
#cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
#cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
#ax.view_init(elev = 17, azim = 100)
#savefig('plots/grad_A0_93_3.png')#, dpi = 300)
#pl.show()
#
#fig = pl.figure()
#ax = fig.gca(projection = '3d')
##ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
##pl.figure()
##contourf(X,Y,h, 1000, cmap = hot())
#surf=ax.plot_surface(1e9*LS,NS,X,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
#CS = contour(1e9*LS,NS,X, colors = 'k', linewidth = 0.5)
#cbar = pl.colorbar(surf)
#cbar.ax.set_ylabel(r'$\mathcal{A}^{(0)}_{[9,3]}\,\,\,[zJ]$', size = 14)
#cbar.add_lines(CS)
##cset = ax.contour(1e9*LS,NS,X,zdir='x',offset = -1,cmap = cm.bone)#winter)
#cset = ax.contour(1e9*LS,NS,X,zdir='z',offset = -0.25,cmap = cm.Blues)
#ax.set_xlabel(r'$\rm{separation}\,\,\,\ell\,\,[nm]$',   size = 18)
#ax.set_ylabel(r'$\rm{n^{th}}\,Matsubra\,term$', size = 18)
#ax.set_zlabel(r'$\nabla \mathcal{A}^{(0)}\,\,[zJ]$',size = 18 ,rotation = 'horizontal' )
#pl.title(r'$\rm{Gradient\, \mathcal{A}^{(0)}\, for \,[9,3]\, in\,water}$',size = 21)
#ax.view_init(elev = 17, azim = 130)
#savefig('plots/grad_A0_93_4.png')#, dpi = 300)
#pl.show()
#
##fig = pl.figure()
##ax = fig.gca(projection = '3d')
###ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
###pl.figure()
###contourf(X,Y,h, 1000, cmap = hot())
##surf=ax.plot_surface(1e9*LS,NS,Y,rstride=1,cstride=1,alpha=0.7,cmap=cm.winter,linewidth = 0.05, antialiased = True, shade = False)
##CS = contour(1e9*LS,NS,Y, colors = 'k', linewidth = 0.5)
##pl.show()
#
#pl.figure()
#pl.contourf(1e9*LS,NS,X, cmap = cm.Blues)
##pl.quiver(1e9*LS[:,::10],NS[:,::10],X[:,::10], units = 'x',pivot='mid',width=0.02)
#pl.quiver(1e9*LS,NS,X)#, units = 'x',pivot='mid')#,width=0.02)
##        np.log(A0[:,::10]), units = 'x',pivot='mid',width=0.02)
#pl.show()
##    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
##    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
##pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
##pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
##pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
##pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
##pl.legend(loc = 'best')
###pl.title(r'65w65 Matsubara terms')
###pl.title(r'90w90 Matsubara terms')
###pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms with static term = 0')
###pl.title(r'290w290 Matsubara terms')
##pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
##pl.xlabel(r'$N$')
###pl.savefig('plots/65_A_vs_n.pdf')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_n0_A_vs_n.eps')
###pl.savefig('plots/93_npk_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
##pl.show()
##
##print 'A0(separation) = ',sum_A0
##print 'A2(separation) = ',sum_A2
##print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
##print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
##
###np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
###
##np.savetxt('data/A0_93_n0_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_n0_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_n0_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_93_npk_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_93_npk_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_93_npk_perpendicular_ret.txt',Ls)
###
###np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
###np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
###np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)
##
##A_py_par  = r'$\mathcal{A}_{\parallel}\sf{[python]}$'
##A0_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
##A2_py_per = r'$\mathcal{A}_{perpend}\sf{[python]}$'
##A_GH_par  = r'$\mathcal{A}_{\parallel}\sf{[ G.H. ]}$'
##
##x_ax = r'$\,\ell\,\,\,\rm{[nm]}$'
##y_ax_par = r'$\mathrm{\mathcal{A}_\parallel(\ell)}\,\,\,\rm{[zJ]}$'
##y_ax_per = r'$\mathrm{\mathcal{A}_\perp (\ell)}\,\,\,\rm{[zJ]}$'
##
##def title(cnt1,cnt2,orientation):
##	return r'$\mathrm{[%s,%s]\,\,Hamaker\,coeff:\, %s \,in\,water,\,retarded}$'%(cnt1,cnt2,orientation)
##
##def svfig(cnt1,cnt2,orientation):
##	return 'plots/140322_%sw%s_HCs_%s_ret.eps'%(cnt1,cnt2,orientation)
##
##pl.figure()
##pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.legend(loc = 'best')
##
###pl.title(title('6','5','Log-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
##pl.title(title('9','3','A(n=0)=0, perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
###pl.savefig(svfig('65','65','perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','n0_perpendicular'))
###pl.savefig(svfig('93','93','npk_perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
##pl.show()
##
##pl.figure()
##pl.loglog(2.*Ls/c,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
###pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
###pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.xlabel(r'$2\ell/c$')
###pl.ylabel(r'(1/$\zeta_0$)2l/c')
##pl.legend(loc = 'best')
##
###pl.title(title('6','5','Log-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
##pl.title(title('9','3','A(n=0)=0, perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
###pl.savefig(svfig('65','65','perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','deltaT_n0_perpendicular'))
###pl.savefig(svfig('93','93','npk_perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
##pl.show()
##
##pl.figure()
##pl.plot(2.*np.sqrt(eiz_w[0])*Ls*zs[0]/c,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
###pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
###pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.xlabel(r'$2\ell\sqrt{\varepsilon(i\zeta_{0})}/c$')
###pl.ylabel(r'(1/$\zeta_0$)2l/c')
##pl.legend(loc = 'best')
##
###pl.title(title('6','5','Log-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','(no first peak eps2_z) perpendicular'))
##pl.title(title('9','3','A(n=0)=0, perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
###pl.savefig(svfig('65','65','perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
##pl.savefig(svfig('93','93','p_n0_perpendicular'))
###pl.savefig(svfig('93','93','npk_perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
##pl.show()
##pl.figure()
##pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'      %(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.legend(loc = 'best')
##pl.minorticks_on()
##pl.ticklabel_format(axis = 'both')
##pl.grid(which = 'both')
##pl.tick_params(which = 'both',labelright = True)
##pl.legend(loc = 'best')
###pl.title(title('6','5','Semi-log, perpendicular'))
###pl.title(title('9','3','(no first peak eps2_z) Semi-log, perpendicular'))
##pl.title(title('9','3','A(n=0)=0,Semilog, perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
###pl.savefig(svfig('65','65','semilog_perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
###pl.savefig(svfig('93','93','npk_semilog_perpendicular'))
##pl.savefig(svfig('93','93','n0_semilog_perpendicular'))
##pl.savefig(svfig('290','290','perpendicular'))
#pl.show()
#
#
#
#
#



