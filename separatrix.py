#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v11.0.3'
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
from math import pi,cos,sin,sqrt,radians,degrees
from setutil import DEBUG_ON,DEBUG_OFF
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as C

def w2phi(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,phi):
    """ longitudinal acceptance: T.Wangler (6.47-48) pp.185 """
    A = 2*pi/(beta*gamma)**3/lamb
    B = q*Ez0*ttf/m0c2
    H = -B*(sin(phisync)-phisync*cos(phisync))
    w2 = 2./A*(H-B*(sin(phi)-phi*cos(phisync)))
    return w2

def test0():
    def sepplot(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,xarr):
        abscisse  = []
        ordinatep = []
        ordinatem = []
        for phi in xarr:
            w2 = w2phi(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,phi)
            if w2 >= 0:
                abscisse.append(degrees(phi))   # abscissa in deg
                DW = sqrt(w2)*m0c2
                ordinatep.append(+DW)         # ordinata in MeV
                ordinatem.append(-DW)
        return (abscisse,ordinatep,ordinatem)
    
    q=1
    phisync=radians(-40.)
    Ez0=1.0
    ttf=0.8
    tkin=6.
    m0c2=C.value('proton mass energy equivalent in MeV')
    gamma=1.+tkin/m0c2
    beta=sqrt(1.-1./gamma**2)
    freq=750.e6
    lamb=C.c/freq

    fig,ax = plt.subplots()
    ax.set_ylabel("\u0394W MeV")
    ax.set_xlabel("\u0394\u03A6 degrees")

    # scatter points
    sigma_phi = radians(12.)     # 1 sigma beam
    sigma_DW  = tkin*0.8/100.    # 1 sigma beam
    x = degrees(sigma_phi) * np.random.randn(1750) + degrees(phisync)  # abscissa in deg
    y = sigma_DW * np.random.randn(1750)                            # ordinate in MeV
    plt.scatter(x,y,s=2)
    str2='\u03C3($\phi$)={:.1f} deg, \u03C3(\u0394W)={:.1f} %'.format(degrees(sigma_phi) ,sigma_DW/tkin*100)

    # 1st separatrix
    start = (phisync-2.*abs(phisync))
    stop =  (phisync+2.*abs(phisync))
    step =  (stop-start)/150
    xarr = [x for x in np.arange(start,stop,step)]

    DEBUG_OFF(f'w2phi {(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,phisync)}')
    abscisse,ordinatep,ordinatem = sepplot(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,xarr)

    DEBUG_ON(f'Wmax={m0c2*sqrt(w2phi(q,m0c2,Ez0,ttf,gamma,beta,lamb,phisync,phisync)):.4e} MeV')

    plt.plot(abscisse,ordinatep,color='blue',label='-40 deg')
    plt.plot(abscisse,ordinatem,color='blue')
    str1=f'Ez0={Ez0} MV/m, tk={tkin} MeV, freq={freq*1e-6} MHz, $\phi_{"s"}$={degrees(phisync)} deg, '

    # two more sepa...
    phisy  = [-35,-30.,-25.,-20.]
    colors = ['red','green','orange','black']
    for i in range(len(phisy)):
        c = colors[i]
        abscisse,ordinatep,ordinatem = sepplot(q,m0c2,Ez0,ttf,gamma,beta,lamb,radians(phisy[i]),xarr)
        x = [y for y in map(lambda x: x-(phisy[i]-degrees(phisync)),abscisse)]  # shift abscissa
        plt.plot(x,ordinatep,color=c,label=f'{phisy[i]} deg')
        plt.plot(x,ordinatem,color=c)

    ax.text(0.01,0.035,str1+str2,transform=ax.transAxes,fontsize=8,verticalalignment='top')
    ax.legend()
    plt.show()
if __name__ == '__main__':
    test0()