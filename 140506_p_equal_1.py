#!/usr/bin/env python
import sys
import numpy as np
from scipy.integrate import romb
import matplotlib.pyplot as pl
# use pyreport -l file.py
from pylab import show

eiz_x = np.loadtxt('data/eiz_x_65.txt') # LDS in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_65.txt') # LDS in parallel direction
eiz_w = np.loadtxt('data/eiz_w.txt') # LDS of water, intervening medium

# Constants
c = 2.99e8               # [m/s]
coeff = 2.411e14         # [rad/s]
Temp = 300.              # [K] 
kbT = Temp * 1.3807e-23  # [J]

ns = np.arange(0.,500.) 
zs = ns * coeff         

Ls = np.arange(1e-9,501e-9,10e-9)  # separation distance between 2 cyclinders

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
p1 = np.zeros(shape = (len(Ls),len(ns)))
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

pMin0 = 0.9
pMax0 = 1.1
# Calculate 1 \geq n \leq 500 terms
for i,L in enumerate(Ls):
    print 'Computing A for separation number %3.0f of %3.0f'%(i, len(Ls))
    for j,n in enumerate(ns):
        p[i,j] = Pn(eiz_w[j],zs[j],L)
p_1 = np.ones((len(Ls),len(ns)))
p_2 =  1./ 2.0*np.ones((len(Ls),len(ns)))
p_3 =  1./ 3.0*np.ones((len(Ls),len(ns)))
p_4 =  1./ 4.0*np.ones((len(Ls),len(ns)))
p_5 =  1./ 5.0*np.ones((len(Ls),len(ns)))
p_6 =  1./ 6.0*np.ones((len(Ls),len(ns)))
p_7 =  1./ 7.0*np.ones((len(Ls),len(ns)))
p_8 =  1./ 8.0*np.ones((len(Ls),len(ns)))
p_9 =  1./ 9.0*np.ones((len(Ls),len(ns)))
p_10 = 1./10.0*np.ones((len(Ls),len(ns)))
p_11 = 1./11.0*np.ones((len(Ls),len(ns)))
p_12 = 1./12.0*np.ones((len(Ls),len(ns)))
p_13 = 1./13.0*np.ones((len(Ls),len(ns)))
p_14 = 1./14.0*np.ones((len(Ls),len(ns)))
p_15 = 1./15.0*np.ones((len(Ls),len(ns)))
p_16 = 1./16.0*np.ones((len(Ls),len(ns)))
p_17 = 1./17.0*np.ones((len(Ls),len(ns)))
p_18 = 1./18.0*np.ones((len(Ls),len(ns)))
p_19 = 1./19.0*np.ones((len(Ls),len(ns)))
p_20 = 1./20.0*np.ones((len(Ls),len(ns)))
p_21 = 1./21.0*np.ones((len(Ls),len(ns)))
p_22 = 1./22.0*np.ones((len(Ls),len(ns)))
p_23 = 1./23.0*np.ones((len(Ls),len(ns)))
p_24 = 1./24.0*np.ones((len(Ls),len(ns)))
p_25 = 1./25.0*np.ones((len(Ls),len(ns)))
p_26 = 1./26.0*np.ones((len(Ls),len(ns)))
p_27 = 1./27.0*np.ones((len(Ls),len(ns)))
p_28 = 1./28.0*np.ones((len(Ls),len(ns)))
p_29 = 1./29.0*np.ones((len(Ls),len(ns)))
p_30 = 1./30.0*np.ones((len(Ls),len(ns)))
p_31 = 1./31.0*np.ones((len(Ls),len(ns)))
p_32 = 1./32.0*np.ones((len(Ls),len(ns)))
p_33 = 1./33.0*np.ones((len(Ls),len(ns)))
p_34 = 1./34.0*np.ones((len(Ls),len(ns)))
p_35 = 1./35.0*np.ones((len(Ls),len(ns)))
p_36 = 1./36.0*np.ones((len(Ls),len(ns)))
p_37 = 1./37.0*np.ones((len(Ls),len(ns)))
p_38 = 1./38.0*np.ones((len(Ls),len(ns)))
p_39 = 2.0*np.ones((len(Ls),len(ns)))
#p_39 = 1./39.0*np.ones((len(Ls),len(ns)))
P1  = np.where(np.isclose(p,p_1,1e-4,1e-1))
P2  = np.where(np.isclose(p,p_2,1e-4,1e-1))
P3  = np.where(np.isclose(p,p_3,1e-4,1e-1))
P4  = np.where(np.isclose(p,p_4,1e-4,1e-1))
P5  = np.where(np.isclose(p,p_5,1e-4,1e-1))
P6  = np.where(np.isclose(p,p_6,1e-4,1e-1))
P7  = np.where(np.isclose(p,p_7,1e-4,1e-1))
P8  = np.where(np.isclose(p,p_8,1e-4,1e-1))
P9  = np.where(np.isclose(p,p_9,1e-4,1e-1))
P10 = np.where(np.isclose(p,p_10 ,1e-4,1e-1))
P11 = np.where(np.isclose(p,p_11 ,1e-4,1e-1))
P12 = np.where(np.isclose(p,p_12 ,1e-4,1e-1))
P13 = np.where(np.isclose(p,p_13 ,1e-4,1e-1))
P14 = np.where(np.isclose(p,p_14 ,1e-4,1e-1))
P15 = np.where(np.isclose(p,p_15 ,1e-4,1e-1))
P16 = np.where(np.isclose(p,p_16 ,1e-4,1e-1))
P17 = np.where(np.isclose(p,p_17 ,1e-4,1e-1))
P18 = np.where(np.isclose(p,p_18 ,1e-4,1e-1))
P19 = np.where(np.isclose(p,p_19 ,1e-4,1e-1))
P20 = np.where(np.isclose(p,p_20 ,1e-4,1e-1))
P21 = np.where(np.isclose(p,p_21 ,1e-4,1e-1))
P22 = np.where(np.isclose(p,p_22 ,1e-4,1e-1))
P23 = np.where(np.isclose(p,p_23 ,1e-4,1e-1))
P24 = np.where(np.isclose(p,p_24 ,1e-4,1e-1))
P25 = np.where(np.isclose(p,p_25 ,1e-4,1e-1))
P26 = np.where(np.isclose(p,p_26 ,1e-4,1e-1))
P27 = np.where(np.isclose(p,p_27 ,1e-4,1e-1))
P28 = np.where(np.isclose(p,p_28 ,1e-4,1e-1))
P29 = np.where(np.isclose(p,p_29 ,1e-4,1e-1))
P30 = np.where(np.isclose(p,p_30 ,1e-4,1e-1))
P31 = np.where(np.isclose(p,p_31 ,1e-4,1e-1))
P32 = np.where(np.isclose(p,p_32 ,1e-4,1e-1))
P33 = np.where(np.isclose(p,p_33 ,1e-4,1e-1))
P34 = np.where(np.isclose(p,p_34 ,1e-4,1e-1))
P35 = np.where(np.isclose(p,p_35 ,1e-4,1e-1))
P36 = np.where(np.isclose(p,p_36 ,1e-4,1e-1))
P37 = np.where(np.isclose(p,p_37 ,1e-4,1e-1))
P38 = np.where(np.isclose(p,p_38 ,1e-4,1e-1))
P39 = np.where(np.isclose(p,p_39 ,1e-4,1e-1))
print 'Ls indices = ', P1[0]
print 'ns indices = ', P1[1]
LL = Ls[P1[0]]
NN = ns[P1[1]]

LL01 = LL[:22]
LL02 = LL[22:32]
LL03 = LL[32:39]
LL04 = LL[39:44]
LL05 = LL[44:48]
LL06 = LL[48:52]
LL07 = LL[52:55]
LL08 = LL[55:57]
LL09 = LL[57:59]
LL10 = LL[59:61]
LL11 = LL[61:63]
LL12 = LL[63:65]
LL13 = LL[65:66]
LL14 = LL[66:68]
LL15 = LL[68:94]

NN01 = NN[:22]
NN02 = NN[22:32]
NN03 = NN[32:39]
NN04 = NN[39:44]
NN05 = NN[44:48]
NN06 = NN[48:52]
NN07 = NN[52:55]
NN08 = NN[55:57]
NN09 = NN[57:59]
NN10 = NN[59:61]
NN11 = NN[61:63]
NN12 = NN[63:65]
NN13 = NN[65:66]
NN14 = NN[66:68]
NN15 = NN[68:94]

pl.figure()
pl.plot(ns[P1[1]],1e9*Ls[P1[0]],  'b x')
pl.plot(ns[P2[1]],1e9*Ls[P2[0]],  'g x')
pl.plot(ns[P3[1]],1e9*Ls[P3[0]],  'r x')
pl.plot(ns[P4[1]],1e9*Ls[P4[0]],  'c x')
pl.plot(ns[P5[1]],1e9*Ls[P5[0]],  'y x')
pl.plot(ns[P6[1]],1e9*Ls[P6[0]],  'm x')
pl.plot(ns[P7[1]],1e9*Ls[P7[0]],  'k x')
pl.plot(ns[P8[1]],1e9*Ls[P8[0]],  'b x')
pl.plot(ns[P9[1]],1e9*Ls[P9[0]],  'g x')
pl.plot(ns[P10[1]],1e9*Ls[P10[0]],'r x')
pl.plot(ns[P11[1]],1e9*Ls[P11[0]],'c x')
pl.plot(ns[P12[1]],1e9*Ls[P12[0]],'y x')
pl.plot(ns[P13[1]],1e9*Ls[P13[0]],'m x')
pl.plot(ns[P14[1]],1e9*Ls[P14[0]],'k x')
pl.plot(ns[P15[1]],1e9*Ls[P15[0]],'b x')
pl.plot(ns[P16[1]],1e9*Ls[P16[0]],'g x')
pl.plot(ns[P17[1]],1e9*Ls[P17[0]],'r x')
pl.plot(ns[P18[1]],1e9*Ls[P18[0]],'c x')
pl.plot(ns[P19[1]],1e9*Ls[P19[0]],'y x')
pl.plot(ns[P20[1]],1e9*Ls[P20[0]],'m x')
pl.plot(ns[P21[1]],1e9*Ls[P21[0]],'k x')
pl.plot(ns[P22[1]],1e9*Ls[P22[0]],'b x')
pl.plot(ns[P23[1]],1e9*Ls[P23[0]],'g x')
pl.plot(ns[P24[1]],1e9*Ls[P24[0]],'r x')
pl.plot(ns[P25[1]],1e9*Ls[P25[0]],'c x')
pl.plot(ns[P26[1]],1e9*Ls[P26[0]],'y x')
pl.plot(ns[P27[1]],1e9*Ls[P27[0]],'m x')
pl.plot(ns[P28[1]],1e9*Ls[P28[0]],'k x')
pl.plot(ns[P29[1]],1e9*Ls[P29[0]],'b x')
pl.plot(ns[P30[1]],1e9*Ls[P30[0]],'m x')
pl.plot(ns[P31[1]],1e9*Ls[P31[0]],'k x')
pl.plot(ns[P32[1]],1e9*Ls[P32[0]],'b x')
pl.plot(ns[P33[1]],1e9*Ls[P33[0]],'g x')
pl.plot(ns[P34[1]],1e9*Ls[P34[0]],'r x')
pl.plot(ns[P35[1]],1e9*Ls[P35[0]],'c x')
pl.plot(ns[P36[1]],1e9*Ls[P36[0]],'y x')
pl.plot(ns[P37[1]],1e9*Ls[P37[0]],'m x')
pl.plot(ns[P38[1]],1e9*Ls[P38[0]],'k x')
pl.plot(ns[P39[1]],1e9*Ls[P39[0]],'b x')
pl.show()

pl.figure()
pl.plot(ns[P1[1]],1e9*Ls[P1[0]],  'b .')
pl.plot(ns[P2[1]],1e9*Ls[P2[0]],  'b .')
pl.plot(ns[P3[1]],1e9*Ls[P3[0]],  'b .')
pl.plot(ns[P4[1]],1e9*Ls[P4[0]],  'b .')
pl.plot(ns[P5[1]],1e9*Ls[P5[0]],  'b .')
pl.plot(ns[P6[1]],1e9*Ls[P6[0]],  'b .')
pl.plot(ns[P7[1]],1e9*Ls[P7[0]],  'b .')
pl.plot(ns[P8[1]],1e9*Ls[P8[0]],  'b .')
pl.plot(ns[P9[1]],1e9*Ls[P9[0]],  'b .')
pl.plot(ns[P10[1]],1e9*Ls[P10[0]],'b .')
pl.plot(ns[P11[1]],1e9*Ls[P11[0]],'b .')
pl.plot(ns[P12[1]],1e9*Ls[P12[0]],'b .')
pl.plot(ns[P13[1]],1e9*Ls[P13[0]],'b .')
pl.plot(ns[P14[1]],1e9*Ls[P14[0]],'b .')
pl.plot(ns[P15[1]],1e9*Ls[P15[0]],'b .')
pl.plot(ns[P16[1]],1e9*Ls[P16[0]],'b .')
pl.plot(ns[P17[1]],1e9*Ls[P17[0]],'b .')
pl.plot(ns[P18[1]],1e9*Ls[P18[0]],'b .')
pl.plot(ns[P19[1]],1e9*Ls[P19[0]],'b .')
pl.plot(ns[P20[1]],1e9*Ls[P20[0]],'b .')
pl.plot(ns[P21[1]],1e9*Ls[P21[0]],'b .')
pl.plot(ns[P22[1]],1e9*Ls[P22[0]],'b .')
pl.plot(ns[P23[1]],1e9*Ls[P23[0]],'b .')
pl.plot(ns[P24[1]],1e9*Ls[P24[0]],'b .')
pl.plot(ns[P25[1]],1e9*Ls[P25[0]],'b .')
pl.plot(ns[P26[1]],1e9*Ls[P26[0]],'b .')
pl.plot(ns[P27[1]],1e9*Ls[P27[0]],'b .')
pl.plot(ns[P28[1]],1e9*Ls[P28[0]],'b .')
pl.plot(ns[P29[1]],1e9*Ls[P29[0]],'b .')
pl.plot(ns[P30[1]],1e9*Ls[P30[0]],'b .')
pl.plot(ns[P31[1]],1e9*Ls[P31[0]],'b .')
pl.plot(ns[P32[1]],1e9*Ls[P32[0]],'b .')
pl.plot(ns[P33[1]],1e9*Ls[P33[0]],'b .')
pl.plot(ns[P34[1]],1e9*Ls[P34[0]],'b .')
pl.plot(ns[P35[1]],1e9*Ls[P35[0]],'b .')
pl.plot(ns[P36[1]],1e9*Ls[P36[0]],'b .')
pl.plot(ns[P37[1]],1e9*Ls[P37[0]],'b .')
pl.plot(ns[P38[1]],1e9*Ls[P38[0]],'b .')
pl.plot(ns[P39[1]],1e9*Ls[P39[0]],'b .')
pl.show()

pl.figure()
pl.plot(NN01,LL01,'b *',  
NN02,LL02,'g s',  
NN03,LL03,'r o',  
NN04,LL04,'c d',  
NN05,LL05,'y +',  
NN06,LL06,'m x',  
NN07,LL07,'k *',  
NN08,LL08,'b s',  
NN09,LL09,'g o',  
NN10,LL10,'r d',  
NN11,LL11,'c +',  
NN12,LL12,'y x',  
NN13,LL13,'m *',  
NN14,LL14,'b s',  
NN15,LL15,'k o') 
pl.show()





#p1= [p[i,j]==1.]
        #p11 = p[p<1.001]#= np.nan #NOTE: remove me later
        #p09 = p11[p11>0.999]#= np.nan #NOTE: remove me later
        #print i,j,p09
        #L11 = Ls[p<1.1]#[p[i,:]==p09[i,:]]#= np.nan #NOTE: remove me later
        #n11 = ns[p[i,j]<1.1]#[p==p09]#= np.nan #NOTE: remove me later
        #p1=p[(p[i,j]>pMin0)*(p[i,j]<pMax0)]
        #NS1=ns[(p[i,j]>pMin0)*(p[i,j]<pMax0)]
        #LS1=Ls[(p[i,j]>pMin0)*(p[i,j]<pMax0)]
        #L09 = L11[p11>0.9]#[p[i,:]==p09[i,:]]#= np.nan #NOTE: remove me later
        #n09 = n11[p11[i,j]>0.9]#[p==p09]#= np.nan #NOTE: remove me later
#        # Integrand A0
#        f0 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))/(T*T+1.)\
#                *(2.*(1.+3.*a[j])*(1.+3.*a[j])*T*T*T*T\
#                + 4.*(1.+2.*a[j]+2.*a[j]+3.*a[j]*a[j])*T*T \
#                + 4.*(1.+a[j])*(1.+a[j]))
#                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
#        # Integrand A2                
#        f2 = T*np.exp(-2.*p[i,j]*np.sqrt(T*T+1.))\
#                /(T*T+1.)\
#                *((T*T*T*T +4.*(T*T)+4.)*(1.-a[j])*(1.-a[j]))
#                #*((eiz_z[j]-eiz_w[j])/(eiz_z[j]-eiz_w[j]))\
#        #Ft0 = np.sum(f0)
#        #Ft2 = np.sum(f2)
#        Ft0 = romb(f0)
#        Ft2 = romb(f2)
#        #Ft = romb(f , axis = 1)
#        #Fty =romb(Ft, axis = 0)
#
#        A0[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft0
#        A0[i,0] = (1./2) * delta[0]*delta[0]*Ft0_term0
#        #A0[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft0_term0
#        
#        A2[i,j] = delta[j]*delta[j]*p[i,j]*p[i,j]*p[i,j]*p[i,j]*Ft2
#        A2[i,0] = (1./2) * delta[0]*delta[0]*Ft2_term0
#        #A2[i,0] = 0.#(1./2) * delta[0]*delta[0]*Ft2_term0
#    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
#    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
#    print 1e9*Ls[i],1e21*(kbT/32)*A0[i,0]
#    print 1e9*Ls[i],1e21*(kbT/32)*A2[i,0]
#    print '-------------------'
#print sum_A0[0] 
#print sum_A2[0]
#print '******************'
#np.savetxt('data/A0_n_65.txt',A0)
#np.savetxt('data/A2_n_65.txt',A2)
#np.savetxt('data/Ls_n_65.txt',Ls)
#np.savetxt('data/A0_65_sum.txt',sum_A0)
#np.savetxt('data/A2_65_sum.txt',sum_A2)
#
#pl.figure()
#pl.plot(Ls,(kbT/32)*A0[:,  0], 'b:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  0]))
#pl.plot(Ls,(kbT/32)*A0[:,  1], 'g:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  1]))
#pl.plot(Ls,(kbT/32)*A0[:,  2], 'r:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  2]))
#pl.plot(Ls,(kbT/32)*A0[:,  3], 'c:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  3]))
#pl.plot(Ls,(kbT/32)*A0[:,  4], 'y:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  4]))
#pl.plot(Ls,(kbT/32)*A0[:,  5], 'm:' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[  5]))
#pl.plot(Ls,(kbT/32)*A0[:, 20], 'b-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 20]))
#pl.plot(Ls,(kbT/32)*A0[:, 40], 'g-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 40]))
#pl.plot(Ls,(kbT/32)*A0[:, 60], 'r-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 60]))
#pl.plot(Ls,(kbT/32)*A0[:, 80], 'c-.', label = r'$A^{(0)}(n=%2.1f)$'%(ns[ 80]))
#pl.plot(Ls,(kbT/32)*A0[:,100], 'b-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[100]))
#pl.plot(Ls,(kbT/32)*A0[:,200], 'g-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[200]))
#pl.plot(Ls,(kbT/32)*A0[:,300], 'r-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[300]))
#pl.plot(Ls,(kbT/32)*A0[:,400], 'c-' , label = r'$A^{(0)}(n=%2.1f)$'%(ns[400]))
#pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
#pl.legend(loc = 'best')      
#pl.title(r'65w65 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{n}\,\, [zJ]$')
#pl.xlabel(r'$\ell$')
#pl.savefig('plots/65_A_n0to400_vs_l.pdf')
###pl.savefig('plots/90_A_vs_n.pdf')
###pl.savefig('plots/91_A_vs_n.pdf')
###pl.savefig('plots/93_A_vs_n.pdf')
###pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()                    

pl.figure()
pl.plot(Ls,2*p[:,  0], 'b:' , label = r'$p(n=%2.1f)$'%(ns[  0]))
pl.plot(Ls,2*p[:,  1], 'g:' , label = r'$p(n=%2.1f)$'%(ns[  1]))
pl.plot(Ls,2*p[:,  2], 'r:' , label = r'$p(n=%2.1f)$'%(ns[  2]))
pl.plot(Ls,2*p[:,  3], 'c:' , label = r'$p(n=%2.1f)$'%(ns[  3]))
pl.plot(Ls,2*p[:,  4], 'y:' , label = r'$p(n=%2.1f)$'%(ns[  4]))
pl.plot(Ls,2*p[:,  5], 'm:' , label = r'$p(n=%2.1f)$'%(ns[  5]))
pl.plot(Ls,2*p[:, 20], 'b-.', label = r'$p(n=%2.1f)$'%(ns[ 20]))
pl.plot(Ls,2*p[:, 40], 'g-.', label = r'$p(n=%2.1f)$'%(ns[ 40]))
pl.plot(Ls,2*p[:, 60], 'r-.', label = r'$p(n=%2.1f)$'%(ns[ 60]))
pl.plot(Ls,2*p[:, 80], 'c-.', label = r'$p(n=%2.1f)$'%(ns[ 80]))
pl.plot(Ls,2*p[:,100], 'b-' , label = r'$p(n=%2.1f)$'%(ns[100]))
pl.plot(Ls,2*p[:,200], 'g-' , label = r'$p(n=%2.1f)$'%(ns[200]))
pl.plot(Ls,2*p[:,300], 'r-' , label = r'$p(n=%2.1f)$'%(ns[300]))
pl.plot(Ls,2*p[:,400], 'c-' , label = r'$p(n=%2.1f)$'%(ns[400]))
#pl.plot(Ls,            sum_A0, 'k-' , label = r'$\Sigma A^{(0)}(n)$'         )            
pl.legend(loc = 'best')      
pl.title(r'65w65 Matsubara terms')
pl.ylabel(r'$\mathcal{A}^{(0)}_{n}\,\, [zJ]$')
pl.xlabel(r'$\ell$')
pl.savefig('plots/65_p_n0to400_vs_l.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
pl.show()                    

pl.figure()
pl.plot(ns,2*p[ 0,:], 'b:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 0]))
pl.plot(ns,2*p[ 1,:], 'g:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 1]))
pl.plot(ns,2*p[ 2,:], 'r:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 2]))
pl.plot(ns,2*p[ 3,:], 'c:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 3]))
pl.plot(ns,2*p[ 4,:], 'y:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 4]))
pl.plot(ns,2*p[ 5,:], 'm:' , label = r'$p(l=%2.1f)$'%(1e9*Ls[ 5]))
pl.plot(ns,2*p[10,:], 'b-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[10]))
pl.plot(ns,2*p[15,:], 'g-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[15]))
pl.plot(ns,2*p[20,:], 'r-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[20]))
pl.plot(ns,2*p[25,:], 'c-.', label = r'$p(l=%2.1f)$'%(1e9*Ls[25]))
pl.plot(ns,2*p[30,:], 'b-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[30]))
pl.plot(ns,2*p[35,:], 'g-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[35]))
pl.plot(ns,2*p[40,:], 'r-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[40]))
pl.plot(ns,2*p[45,:], 'c-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[45]))
pl.plot(ns,2*p[50,:], 'y-' , label = r'$p(l=%2.1f)$'%(1e9*Ls[50]))
pl.plot(ns,p[10,:]/p[10,:], 'k-')# , label = r'$p(l=%2.1f)$'%(Ls[45]))
pl.legend(loc = 'best')                               
pl.axis([0,500,0,1.2])
pl.title(r'65w65 p(n,l)')
pl.ylabel(r'$p(n,l)$')
pl.xlabel(r'$n$')
pl.savefig('plots/65_p_l_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
pl.show()                    

##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_290_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_290_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_290_perpendicular_ret.txt',Ls)

#        #A0[A0>1e6]= np.nan #NOTE: remove me later
#        #A2[A2>1e6]= np.nan #NOTE: remove me later
#    sum_A0 = (kbT/(32.)) * np.sum(A0, axis = 1)
#    sum_A2 = (kbT/(32.)) * np.sum(A2, axis = 1)
#pl.figure()
#pl.loglog(ns,(kbT/(32.)) * A0[0,:], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A0[1,:], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A0[2,:], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A0[3,:], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.loglog(ns,(kbT/(32.)) * A2[0,:], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(ns,(kbT/(32.)) * A2[1,:], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(ns,(kbT/(32.)) * A2[2,:], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(ns,(kbT/(32.)) * A2[3,:], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
##pl.title(r'290w290 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
#pl.xlabel(r'$N$')
#pl.savefig('plots/65_A_vs_n.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()
#
#pl.figure()
#pl.loglog(p[0,:],(kbT/(32.)) * A0[0,:] / A0[0,0], 'b-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(p[1,:],(kbT/(32.)) * A0[1,:] / A0[1,0], 'g-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(p[2,:],(kbT/(32.)) * A0[2,:] / A0[2,0], 'r-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(p[3,:],(kbT/(32.)) * A0[3,:] / A0[3,0], 'y-', label = r'$A^{0,2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.loglog(p[0,:],(kbT/(32.)) * A2[0,:] / A2[0,0], 'b:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[0]))
#pl.loglog(p[1,:],(kbT/(32.)) * A2[1,:] / A2[1,0], 'g:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[1]))
#pl.loglog(p[2,:],(kbT/(32.)) * A2[2,:] / A2[2,0], 'r:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[2]))
#pl.loglog(p[3,:],(kbT/(32.)) * A2[3,:] / A2[3,0], 'y:')#, label = r'$A^{2}(\ell=%2.1f\,nm)$'%(1e9*Ls[3]))
#pl.legend(loc = 'best')
#pl.title(r'65w65 Matsubara terms vs p')
##pl.title(r'90w90 Matsubara terms')
##pl.title(r'91w91 Matsubara terms')
##pl.title(r'93w93 Matsubara terms')
##pl.title(r'290w290 Matsubara terms')
#pl.ylabel(r'$\mathcal{A}^{(0)}_{N}, \,\, \mathcal{A}^{(2)}_{N}$')
#pl.xlabel(r'$N$')
#pl.savefig('plots/65_A_vs_p.pdf')
##pl.savefig('plots/90_A_vs_n.pdf')
##pl.savefig('plots/91_A_vs_n.pdf')
##pl.savefig('plots/93_A_vs_n.pdf')
##pl.savefig('plots/290_A_vs_n.pdf')
#pl.show()
#
#print 'A0(separation) = ',sum_A0
#print 'A2(separation) = ',sum_A2
#print 'Contribution to A0 from n=0 term = ', (kbT/(12.*np.pi))*A0[:,0]
#print 'Contribution to A2 from n=0 term = ', (kbT/(12.*np.pi))*A2[:,0]
#
#np.savetxt('data/A0_65_perpendicular_ret.txt',sum_A0)
#np.savetxt('data/A2_65_perpendicular_ret.txt',sum_A2)
#np.savetxt('data/Lengths_65_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_90_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_90_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_90_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_91_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_91_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_91_perpendicular_ret.txt',Ls)
##
##np.savetxt('data/A0_93_perpendicular_ret.txt',sum_A0)
##np.savetxt('data/A2_93_perpendicular_ret.txt',sum_A2)
##np.savetxt('data/Lengths_93_perpendicular_ret.txt',Ls)
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
#	return 'plots/140322_%sw%s_HCs_%s_ret.pdf'%(cnt1,cnt2,orientation)
#
##pl.figure()
##pl.loglog(1e9*Ls,1e21*sum_A0,      label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.loglog(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1f nm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.legend(loc = 'best')
##
##pl.title(title('6','5','Log-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
##pl.savefig(svfig('65','65','perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
###pl.savefig(svfig('93','93','perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
##pl.show()
##
##pl.figure()
##pl.semilogy(1e9*Ls,1e21*sum_A0,label=r'$\mathcal{A^{(0)}}(\ell=%1.1f nm)=%3.2f,\,\,\, \mathcal{A^{(0)}}(\ell=%1.1fnm)=%3.2f$'%(1e9*Ls[0],      1e21*sum_A0[0],1e9*Ls[1],1e21*sum_A0[1]))
##pl.semilogy(1e9*Ls,1e21*sum_A2,label=r'$\mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f, \,\,\, \mathcal{A^{(2)}}(\ell=\, %1.1fnm)=%3.2f$'%(1e9*Ls[0],1e21*sum_A2[0],1e9*Ls[1],1e21*sum_A2[1]))
##pl.xlabel(x_ax)
##pl.ylabel(y_ax_per)
##pl.title(title('6','5','perpendicular'))
##pl.legend(loc = 'best')
##pl.minorticks_on()
##pl.ticklabel_format(axis = 'both')
##pl.grid(which = 'both')
##pl.tick_params(which = 'both',labelright = True)
##pl.legend(loc = 'best')
##pl.title(title('6','5','Semi-log, perpendicular'))
###pl.title(title('9','0','perpendicular'))
###pl.title(title('9','1','perpendicular'))
###pl.title(title('9','3','perpendicular'))
###pl.title(title('29','0','perpendicular'))
###
##pl.savefig(svfig('65','65','semilog_perpendicular'))
###pl.savefig(svfig('90','90','perpendicular'))
###pl.savefig(svfig('91','91','perpendicular'))
###pl.savefig(svfig('93','93','perpendicular'))
###pl.savefig(svfig('290','290','perpendicular'))
##pl.show()
##
#

