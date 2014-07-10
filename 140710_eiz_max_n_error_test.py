#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show
import pyreport

#eV_x,eps2_x = np.loadtxt('data/CNT6_5_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT6_5_ze2_solid_30.txt',unpack=True, usecols = [0,1])

#eV_x,eps2_x = np.loadtxt('data/CNT9_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT9_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#eps2_z[:57] = 0.0 # removes first divergent peak in eps2_z, saves file with npk

#eV_x,eps2_x = np.loadtxt('data/CNT9_1_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT9_1_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#
#eV_x,eps2_x = np.loadtxt('data/CNT9_3_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT9_3_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#eps2_z[:120] = 0.0 # removes first divergent peak in eps2_z, saves file with npk
#
#eV_x,eps2_x = np.loadtxt('data/CNT29_0_xe2_solid_30.txt',unpack=True, usecols = [0,1])
#eV_z,eps2_z = np.loadtxt('data/CNT29_0_ze2_solid_30.txt',unpack=True, usecols = [0,1])

eV_x,eps2_x = np.loadtxt('data/CNT24_24_xe2_solid_30.txt',unpack=True, usecols = [0,1])
eV_z,eps2_z = np.loadtxt('data/CNT24_24_ze2_solid_30.txt',unpack=True, usecols = [0,1])
#
#eV_w,eps2_w = np.loadtxt('data/water-L.txt',unpack=True, usecols = [0,1])

#$ Make a list of Matsubara energies as integer multiples of N
#$ \xi_{N}^(RT) = (2 \pi k_{B} T / hbar)N = coeff*N

eV_w, eps2_w  = np.loadtxt('data/water_fine_DD_140708.txt', unpack=True, usecols = [0,1])
#x    ,        y = np.loadtxt('data/water_finer_DD_140708.txt', unpack=True,usecols =[0,1]) # this is the LDS

coeff = 0.159      # [eV] 
T = 300.0          # Room temp
zs20 = coeff*np.arange(1,2000)      # energies 
zs10 = coeff*np.arange(1,1000)      # energies 
zs05 = coeff*np.arange(1,500)      # energies 
zs02 = coeff*np.arange(1,200)      # energies 

#$ Make empty lists for \varepsilon(i\xi) and integrand
eiz_x20 = np.empty(len(zs20))
eiz_z20 = np.empty(len(zs20))
eiz_w20 = np.empty(len(zs20))

eiz_x10 = np.empty(len(zs10))
eiz_z10 = np.empty(len(zs10))
eiz_w10 = np.empty(len(zs10))

eiz_x05 = np.empty(len(zs05))
eiz_z05 = np.empty(len(zs05))
eiz_w05 = np.empty(len(zs05))

eiz_x02 = np.empty(len(zs02))
eiz_z02 = np.empty(len(zs02))
eiz_w02 = np.empty(len(zs02))

integrand_x20=np.empty(len(eV_x))
integrand_z20=np.empty(len(eV_z))
integrand_w20=np.empty(len(eV_w))

integrand_x10=np.empty(len(eV_x))
integrand_z10=np.empty(len(eV_z))
integrand_w10=np.empty(len(eV_w))

integrand_x05=np.empty(len(eV_x))
integrand_z05=np.empty(len(eV_z))
integrand_w05=np.empty(len(eV_w))

integrand_x02=np.empty(len(eV_x))
integrand_z02=np.empty(len(eV_z))
integrand_w02=np.empty(len(eV_w))

#$ Compute \varepsilon(i\xi) for each \xi_N
for j in range(len(zs20)):
    for k in range(len(eV_x)):
        integrand_x20[k]=eV_x[k]*eps2_x[k] / (eV_x[k]**2 + zs20[j]**2)
    eiz_x20[j] = 1 + (2./np.pi) * np.trapz(integrand_x20,eV_x)

    for m in range(len(eV_z)):
        integrand_z20[m]=eV_z[m]*eps2_z[m] / (eV_z[m]**2 + zs20[j]**2)
    eiz_z20[j] = 1 + (2./np.pi) * np.trapz(integrand_z20,eV_z)    

    for p in range(len(eV_w)):
        integrand_w02[p]=eV_w[p]*eps2_w[p] / (eV_w[p]**2 + zs20[j]**2)
    eiz_w20[j] = 1 + (2./np.pi) * np.trapz(integrand_w20,eV_w)    

for j in range(len(zs10)):
    for k in range(len(eV_x)):
        integrand_x10[k]=eV_x[k]*eps2_x[k] / (eV_x[k]**2 + zs10[j]**2)
    eiz_x10[j] = 1 + (2./np.pi) * np.trapz(integrand_x10,eV_x)

    for m in range(len(eV_z)):
        integrand_z10[m]=eV_z[m]*eps2_z[m] / (eV_z[m]**2 + zs10[j]**2)
    eiz_z10[j] = 1 + (2./np.pi) * np.trapz(integrand_z10,eV_z)    

    for p in range(len(eV_w)):
        integrand_w10[p]=eV_w[p]*eps2_w[p] / (eV_w[p]**2 + zs10[j]**2)
    eiz_w10[j] = 1 + (2./np.pi) * np.trapz(integrand_w10,eV_w)    

for j in range(len(zs05)):
    for k in range(len(eV_x)):
        integrand_x05[k]=eV_x[k]*eps2_x[k] / (eV_x[k]**2 + zs05[j]**2)
    eiz_x05[j] = 1 + (2./np.pi) * np.trapz(integrand_x05,eV_x)

    for m in range(len(eV_z)):
        integrand_z05[m]=eV_z[m]*eps2_z[m] / (eV_z[m]**2 + zs05[j]**2)
    eiz_z05[j] = 1 + (2./np.pi) * np.trapz(integrand_z05,eV_z)    

    for p in range(len(eV_w)):
        integrand_w05[p]=eV_w[p]*eps2_w[p] / (eV_w[p]**2 + zs05[j]**2)
    eiz_w05[j] = 1 + (2./np.pi) * np.trapz(integrand_w05,eV_w)    

for j in range(len(zs02)):
    for k in range(len(eV_x)):
        integrand_x02[k]=eV_x[k]*eps2_x[k] / (eV_x[k]**2 + zs02[j]**2)
    eiz_x02[j] = 1 + (2./np.pi) * np.trapz(integrand_x02,eV_x)

    for m in range(len(eV_z)):
        integrand_z02[m]=eV_z[m]*eps2_z[m] / (eV_z[m]**2 + zs02[j]**2)
    eiz_z02[j] = 1 + (2./np.pi) * np.trapz(integrand_z02,eV_z)    

    for p in range(len(eV_w)):
        integrand_w02[p]=eV_w[p]*eps2_w[p] / (eV_w[p]**2 + zs02[j]**2)
    eiz_w02[j] = 1 + (2./np.pi) * np.trapz(integrand_w02,eV_w)    

def Delta_par(par,med):
	return ((par-1.) - (med-1.))/(med-1.)

def Delta_prp(perp,med):
	return ((perp-1.) - (med-1.))/((perp-1.) + (med-1.))


dpar20 = Delta_par(eiz_z20/eiz_z20[0],eiz_w20/eiz_w20[0]) 
dpar10 = Delta_par(eiz_z10/eiz_z10[0],eiz_w10/eiz_w10[0]) 
dpar05 = Delta_par(eiz_z05/eiz_z05[0],eiz_w05/eiz_w05[0]) 
dpar02 = Delta_par(eiz_z02/eiz_z02[0],eiz_w02/eiz_w02[0]) 

dprp20 = Delta_prp(eiz_x20/eiz_x20[0],eiz_w20/eiz_w20[0]) 
dprp10 = Delta_prp(eiz_x10/eiz_x10[0],eiz_w10/eiz_w10[0]) 
dprp05 = Delta_prp(eiz_x05/eiz_x05[0],eiz_w05/eiz_w05[0]) 
dprp02 = Delta_prp(eiz_x02/eiz_x02[0],eiz_w02/eiz_w02[0]) 

Dpar20 = np.trapz(dpar20)# = Delta_par(eiz_z20,eiz_w20) 
Dpar10 = np.trapz(dpar10)# = Delta_par(eiz_z10,eiz_w10) 
Dpar05 = np.trapz(dpar05)# = Delta_par(eiz_z05,eiz_w05) 
Dpar02 = np.trapz(dpar02)# = Delta_par(eiz_z02,eiz_w02) 

Dprp20 = np.trapz(dprp20)# = Delta_prp(eiz_x20,eiz_w20) 
Dprp10 = np.trapz(dprp10)# = Delta_prp(eiz_x10,eiz_w10) 
Dprp05 = np.trapz(dprp05)# = Delta_prp(eiz_x05,eiz_w05) 
Dprp02 = np.trapz(dprp02)# = Delta_prp(eiz_x02,eiz_w02) 

print 'Dpar20 = ',Dpar20
print 'Dpar10 = ',Dpar10
print 'Dpar05 = ',Dpar05
print 'Dpar02 = ',Dpar02

print 'Dprp20 = ',Dprp20
print 'Dprp10 = ',Dprp10
print 'Dprp05 = ',Dprp05
print 'Dprp02 = ',Dprp02






pl.figure()
pl.plot(zs20, dpar20,'b-', label = r'$\Delta_{\parallel}(n_{max}=2,000)$')
pl.plot(zs20, dprp20,'b--',label = r'$\Delta_{\perp}(n_{max}=2,000)$')
pl.plot(zs10, dpar10,'g-', label = r'$\Delta_{\parallel}(n_{max}=1,000)$')
pl.plot(zs10, dprp10,'g--',label = r'$\Delta_{\perp}(n_{max}=1,000)$')
pl.plot(zs05, dpar05,'r-', label = r'$\Delta_{\parallel}(n_{max}=500)$')
pl.plot(zs05, dprp05,'r--',label = r'$\Delta_{\perp}(n_{max}=500)$')
pl.plot(zs02, dpar02,'c-', label = r'$\Delta_{\parallel}(n_{max}=200)$')
pl.plot(zs02, dprp02,'c--',label = r'$\Delta_{\perp}(n_{max}=200)$')
pl.plot(loc = 'best')
pl.xlabel(r'$\omega_n$', size = 24)
pl.ylabel(r'$\epsilon(i\omega_{n})$' , size = 24)
pl.savefig('plots/24_error_due_to_n.pdf')
pl.show()

pl.figure()
pl.plot(np.max(zs20),100.* np.trapz(eiz_x20)/np.trapz(eiz_x20),'b *', label = r'$\varepsilon_{\hat{x}}(n_{max}=2,000)$')
pl.plot(np.max(zs20),100.* np.trapz(eiz_z20)/np.trapz(eiz_z20),'r *', label = r'$\varepsilon_{\hat{z}}(n_{max}=2,000)$')
pl.plot(np.max(zs20),100.* np.trapz(eiz_w20)/np.trapz(eiz_w20),'c *', label = r'$\varepsilon_{\hat{w}}(n_{max}=2,000)$')

pl.plot(np.max(zs10),100.*(np.trapz(eiz_x20)-np.trapz(eiz_x10))/np.trapz(eiz_x20),'b +', label = r'$\varepsilon_{\hat{x}}(n_{max}=1,000)$')
pl.plot(np.max(zs10),100.*(np.trapz(eiz_z20)-np.trapz(eiz_z10))/np.trapz(eiz_z20),'r +', label = r'$\varepsilon_{\hat{z}}(n_{max}=1,000)$')
pl.plot(np.max(zs10),100.*(np.trapz(eiz_w20)-np.trapz(eiz_w10))/np.trapz(eiz_w20),'c +', label = r'$\varepsilon_{\hat{w}}(n_{max}=1,000)$')

pl.plot(np.max(zs05),100.*(np.trapz(eiz_x20)-np.trapz(eiz_x05))/np.trapz(eiz_x20),'b x', label = r'$\varepsilon_{\hat{x}}(n_{max}=500)$')
pl.plot(np.max(zs05),100.*(np.trapz(eiz_z20)-np.trapz(eiz_z05))/np.trapz(eiz_z20),'r x', label = r'$\varepsilon_{\hat{z}}(n_{max}=500)$')
pl.plot(np.max(zs05),100.*(np.trapz(eiz_w20)-np.trapz(eiz_w05))/np.trapz(eiz_w20),'c x', label = r'$\varepsilon_{\hat{w}}(n_{max}=500)$')

pl.plot(np.max(zs02),100.*(np.trapz(eiz_x20)-np.trapz(eiz_x02))/np.trapz(eiz_x20),'b ^', label = r'$\varepsilon_{\hat{x}}(n_{max}=200)$')
pl.plot(np.max(zs02),100.*(np.trapz(eiz_z20)-np.trapz(eiz_z02))/np.trapz(eiz_z20),'r ^', label = r'$\varepsilon_{\hat{z}}(n_{max}=200)$')
pl.plot(np.max(zs02),100.*(np.trapz(eiz_w20)-np.trapz(eiz_w02))/np.trapz(eiz_w20),'c ^', label = r'$\varepsilon_{\hat{w}}(n_{max}=200)$')
pl.legend(loc = 'best')
pl.xlabel(r'$\omega_n$', size = 24)
pl.ylabel(r'$\epsilon(i\omega_{n})$' , size = 24)
pl.savefig('plots/24_error_due_to_n.pdf')
pl.show()

pl.figure()
pl.loglog(zs20, eiz_x20,'b-', label = r'$\varepsilon_{\hat{x}}(n_{max}=2,000)$')
pl.loglog(zs20, eiz_z20,'r-', label = r'$\varepsilon_{\hat{z}}(n_{max}=2,000)$')
pl.loglog(zs20, eiz_w20,'c-', label = r'$\varepsilon_{\hat{w}}(n_{max}=2,000)$')

pl.loglog(zs10,(eiz_x10),'b--', label = r'$\varepsilon_{\hat{x}}(n_{max}=1,000)$')
pl.loglog(zs10,(eiz_z10),'r--', label = r'$\varepsilon_{\hat{z}}(n_{max}=1,000)$')
pl.loglog(zs10,(eiz_w10),'c--', label = r'$\varepsilon_{\hat{w}}(n_{max}=1,000)$')

pl.loglog(zs05,(eiz_x05),'b-.', label = r'$\varepsilon_{\hat{x}}(n_{max}=500)$')
pl.loglog(zs05,(eiz_z05),'r-.', label = r'$\varepsilon_{\hat{z}}(n_{max}=500)$')
pl.loglog(zs05,(eiz_w05),'c-.', label = r'$\varepsilon_{\hat{w}}(n_{max}=500)$')

pl.loglog(zs02,(eiz_x02),'b:', label = r'$\varepsilon_{\hat{x}}(n_{max}=200)$')
pl.loglog(zs02,(eiz_z02),'r:', label = r'$\varepsilon_{\hat{z}}(n_{max}=200)$')
pl.loglog(zs02,(eiz_w02),'c:', label = r'$\varepsilon_{\hat{w}}(n_{max}=200)$')
pl.legend(loc = 'best')
pl.xlabel(r'$\omega_n$', size = 24)
pl.ylabel(r'$\epsilon(i\omega_{n})$' , size = 24)
pl.savefig('plots/24_error_due_to_n.pdf')
pl.show()


pl.figure()
pl.plot(np.max(zs20), np.trapz(eiz_x20)/np.trapz(eiz_x20),'b *', label = r'$\varepsilon_{\hat{x}}(n_{max}=2,000)$')
pl.plot(np.max(zs20), np.trapz(eiz_z20)/np.trapz(eiz_z20),'r *', label = r'$\varepsilon_{\hat{z}}(n_{max}=2,000)$')
pl.plot(np.max(zs20), np.trapz(eiz_w20)/np.trapz(eiz_w20),'c *', label = r'$\varepsilon_{\hat{w}}(n_{max}=2,000)$')

pl.plot(np.max(zs10),(np.trapz(eiz_x10))/np.trapz(eiz_x20),'b +', label = r'$\varepsilon_{\hat{x}}(n_{max}=1,000)$')
pl.plot(np.max(zs10),(np.trapz(eiz_z10))/np.trapz(eiz_z20),'r +', label = r'$\varepsilon_{\hat{z}}(n_{max}=1,000)$')
pl.plot(np.max(zs10),(np.trapz(eiz_w10))/np.trapz(eiz_w20),'c +', label = r'$\varepsilon_{\hat{w}}(n_{max}=1,000)$')

pl.plot(np.max(zs05),(np.trapz(eiz_x05))/np.trapz(eiz_x20),'b x', label = r'$\varepsilon_{\hat{x}}(n_{max}=500)$')
pl.plot(np.max(zs05),(np.trapz(eiz_z05))/np.trapz(eiz_z20),'r x', label = r'$\varepsilon_{\hat{z}}(n_{max}=500)$')
pl.plot(np.max(zs05),(np.trapz(eiz_w05))/np.trapz(eiz_w20),'c x', label = r'$\varepsilon_{\hat{w}}(n_{max}=500)$')

pl.plot(np.max(zs02),(np.trapz(eiz_x02))/np.trapz(eiz_x20),'b ^', label = r'$\varepsilon_{\hat{x}}(n_{max}=200)$')
pl.plot(np.max(zs02),(np.trapz(eiz_z02))/np.trapz(eiz_z20),'r ^', label = r'$\varepsilon_{\hat{z}}(n_{max}=200)$')
pl.plot(np.max(zs02),(np.trapz(eiz_w02))/np.trapz(eiz_w20),'c ^', label = r'$\varepsilon_{\hat{w}}(n_{max}=200)$')
pl.legend(loc = 'best')
pl.xlabel(r'$\omega_n$', size = 24)
pl.ylabel(r'$\epsilon(i\omega_{n})$' , size = 24)
pl.savefig('plots/24_error_due_to_n.pdf')
pl.show()
