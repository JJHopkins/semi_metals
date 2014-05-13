In [149]: pl.loglog(Ls,A93l, label = r'$\mathcal{A}^{(0)}(\ell)$')
Out[149]: [<matplotlib.lines.Line2D at 0x1878c4d0>]

In [150]: pl.loglog(Ls,-X93l, label = r'$-d\mathcal{A}/d\ell$')
Out[150]: [<matplotlib.lines.Line2D at 0x199f9c50>]

In [151]: pl.legend(loc = 'best')
Out[151]: <matplotlib.legend.Legend at 0x11a5ced0>

In [152]: pl.xlabel(r'separation $\ell$ [m]')
Out[152]: <matplotlib.text.Text at 0x187656d0>

In [153]: pl.ylabel(r'$\mathcal{A}^{(0)}(\ell),\,\,-d\mathcal{A}/d\ell$')
Out[153]: <matplotlib.text.Text at 0x1964e890>

In [154]: pl.title('Separation value at knee and at max derivative ')
Out[154]: <matplotlib.text.Text at 0x2688ad10>

In [155]: pl.show()

###### delta plots####
NS = np.arange(1,500)
D = np.diff(delta)
D2 = np.diff(delta*delta)
Hline = D/D
Hline0 = D-D

pl.plot(NS,Hline,'k--')
pl.plot(NS,Hline0,'k--')

pl.plot(ns, delta, 'b-', label=r'$\Delta$')
pl.plot(ns, delta*delta, 'g-', label=r'$\Delta^2$')

pl.plot(NS,D,'r-', label = r'$d\Delta/dn$')
pl.plot(NS,D2,'c-', label = r'$d\Delta^2/dn$')

pl.legend(loc = 'best')
pl.xlabel('n')
#pl.ylabel(r'$d\,/dn$')
pl.title('[6,5] Derivative of dielectric contrast with respect to n')
#pl.title('[9,3] Derivative of dielectric contrast with respect to n')
pl.show()


###### delta plots####
NS = np.arange(1,500)
A = np.diff(a)
A2 = np.diff(a*a)
Hline = A/A
Hline0 = A-A

pl.plot(NS,Hline,'k--')
pl.plot(NS,Hline0,'k--')

pl.plot(ns, a, 'b-', label=r'$a(i\zeta_n)$')
pl.plot(ns, a*a, 'g-', label=r'$a(i\zeta_n)^2$')

pl.plot(NS,A,'r-', label = r'$da/dn$')
pl.plot(NS,A2,'c-', label = r'$da^2/dn$')

pl.legend(loc = 'best')
pl.xlabel('n')
#pl.ylabel(r'$d\,/dn$')
#pl.title('[6,5] Derivative of anisotropy ratio with respect to n')
pl.title('[9,3] Derivative of anisotropy ratio with respect to n')
pl.show()


