#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show
import pyreport

#$  Jaime Hopkins HW2 prob 2
#$  % Python code for Kramers-Kronig transform of
#$ \varepsilon\prime\prime(i\omega)

#$  \begin{equation}
#$  p(\cos\theta | b/a) = N ( 1 + (b/a) \cos^2 \theta )
#$  \end{equation}
#$  %
#$  \subsection{First, find N such that p = 1:}
#$  \begin{equation}
#$  p = &&\int \limits_{-1}^{1} N ( 1 + (b/a) \cos^2 \theta ) d (\cos \theta) = 1 
#$  \end{equation}
#$  %
#$  \begin{eqnarray}
#$  p = &&\int \limits_{-1}^{1} N ( 1 + (b/a) \cos^2 \theta ) d (\cos \theta) = 1 
#$  1/N = &&\int \limits_{-1}^{1} ( 1 + (b/a) \cos^2 \theta ) d (\cos \theta)
#$   = &&\int \limits_{-1}^{1} ( 1 + (b/a) \x^2 ) dx = 2(1 + (b/a)/3)
#$  %
#$  N = \frac{3}{2(3+b/a)}
#$  \end{eqnarray}
#$  %
#$ \subsection{2 (a) using the method of moments} 
#$ \[\begin{eqnarray}
#$   \left<\cos^2\theta\right> &= N \int \limits_{-1}^{1} \mathrm{d}\cos\theta\left(\cos^2 \theta + (b/a) \cos^4 \theta\right)\\
#$                             &= N \left.\left(\frac{1}{3}x^3 + \frac{(b/a)}{5}x^5\right)\right|_{-1}^{1}\\
#$                             &= \frac{1}{5}\cdot\frac{5 + 3(b/a)}{3 +(b/a)}
#$  \end{array}
#$  %

#$ First, input \varepsilon\prime\prime(i\omega) files from Dan:

#eV_x,eps2_x = np.loadtxt('data/CNT6_5_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT6_5_ze2_solid_30.txt',unpack=True, usecols = [0,1])

#eV_x,eps2_x = np.loadtxt('data/CNT9_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT9_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#eps2_z[:57] = 0.0 # removes first divergent peak in eps2_z, saves file with npk

#eV_x,eps2_x = np.loadtxt('data/CNT9_1_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT9_1_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#
eV_x,eps2_x = np.loadtxt('data/CNT9_3_xe2_solid_30.txt',unpack=True, usecols = [0,1])
eV_z,eps2_z = np.loadtxt('data/CNT9_3_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#eps2_z[:120] = 0.0 # removes first divergent peak in eps2_z, saves file with npk
#
#eV_x,eps2_x = np.loadtxt('data/CNT29_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT29_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])

eV_w,eps2_w = np.loadtxt('data/water-L.txt',unpack=True, usecols = [0,1])

#$ Make a list of Matsubara energies as integer multiples of N
#$ \xi_{N}^(RT) = (2 \pi k_{B} T / hbar)N = coeff*N
coeff = 0.159      # [eV] 
T = 300.0          # Room temp
n = np.arange(0,500) 
z = n * coeff      # energies 

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

#$ Make empty lists for \varepsilon(i\xi) and integrand
eiz_x = np.empty(len(z))
eiz_z = np.empty(len(z))
eiz_w = np.empty(len(z))

integrand_x=np.empty(len(eV_x))
integrand_z=np.empty(len(eV_z))
integrand_w=np.empty(len(eV_w))

#$ Compute \varepsilon(i\xi) for each \xi_N
for j in range(len(z)):
    for k in range(len(eV_x)):
        integrand_x[k]=eV_x[k]*eps2_x[k] / (eV_x[k]**2 + z[j]**2)
    eiz_x[j] = 1 + (2./np.pi) * np.trapz(integrand_x,eV_x)

    for m in range(len(eV_z)):
        integrand_z[m]=eV_z[m]*eps2_z[m] / (eV_z[m]**2 + z[j]**2)
    eiz_z[j] = 1 + (2./np.pi) * np.trapz(integrand_z,eV_z)    

    for p in range(len(eV_w)):
        integrand_w[p]=eV_w[p]*eps2_w[p] / (eV_w[p]**2 + z[j]**2)
    eiz_w[j] = 1 + (2./np.pi) * np.trapz(integrand_w,eV_w)    

np.savetxt( "data/energies.txt", z    )
np.savetxt( "data/eiz_x_93.txt", eiz_x)
#np.savetxt( "data/eiz_z_npk_93.txt", eiz_z)
np.savetxt( "data/eiz_z_93.txt", eiz_z)
np.savetxt( "data/eiz_w.txt"   , eiz_w)
#
a =  Aiz(eiz_x,eiz_z,eiz_w)
#
fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(n,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
ax.plot(n,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{n})$')
ax.plot(n,eiz_z, color = 'r', label = r'$max\,%6.2f$'%max(eiz_z))
ax.plot(n,eiz_w, color = 'c', label = r'$\varepsilon_{\hat{w}}(i\zeta_{n})$')
pl.axis([0,500,0,10])
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 24)
#pl.legend()
pl.title(r'[9,3] eiz and water eiz')
#pl.title(r'[9,3] eiz (no first peak in eps2_z) and water eiz')

ax_inset = fig.add_axes([0.53,0.50,0.36,0.36])
ax_inset.plot(n, a,'k-.', linewidth = 2)#,label=r'$a(i\xi_{N})$')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$N$', size = 14)
pl.ylabel(r'$a(i\xi)$', size = 14)
pl.savefig('plots/93w93_eiz.eps')
#pl.savefig('plots/93w93_eiz_npk.eps')
pl.show()


pl.figure()
pl.plot(eV_x,eps2_x, color = 'b', label = r'$\varepsilon^{\prime\prime}_\hat{x}(\omega)$')
pl.plot(eV_z,eps2_z, color = 'r', label = r'$\varepsilon^{\prime\prime}_\hat{z}(\omega)$')
pl.plot(eV_z,eps2_z, color = 'r', label = r'$max\,%.8s$'%max(eps2_z))
pl.plot(eV_w,eps2_w, color = 'c', label = r'$\varepsilon^{\prime\prime}_{H_{2}O}(\omega)$')
pl.xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
pl.ylabel(r'$\varepsilon^{\prime\prime}(\omega)$', size = 24)
pl.axis([0,40,0,30])
pl.legend()
pl.title(r'[9,3] eps2 and water eps2')
#pl.title(r'[9,3] eps2 (no first peak) and water eps2')
pl.savefig('plots/93w93_eps2.eps')
#pl.savefig('plots/93w93_npk_eps2.eps')
show()


