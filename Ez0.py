#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
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
import matplotlib.pyplot as plt
import numpy as np
from math import sin,tan,pi,exp,fmod
from collections import namedtuple

from setutil import PARAMS,DEBUG,Proton

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

'''
Cavity E(z,r=0) field profile
'''

def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res

def Kpoly(z,sigma,mu,E):
    """
    Calculate polynom coefficients from NG formula
    """
    poly = []
    Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0','coeff'])
    anz = len(z)
    for i in range(0,anz-2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  = z0-zl
        El  = E*NGauss(zl,sigma,mu)
        E0  = E*NGauss(z0,sigma,mu)
        Er  = E*NGauss(zr,sigma,mu)
        b = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
        a = (Er-El)/(2*E0*dz)           # getetstet mit Bleistift u. Papier
        pval = Polyval(zl,z0,zr,dz,b,a,E0,0.)
        poly.append(pval)

        # USING np.polyfit() liefert gleiches Resultat wie Langrange 3 Punkt
        # x   = np.array((zl,z0,zr))
        # y   = np.array((El,E0,Er))
        # coeff = np.polyfit(x,y,2)
        # b   = coeff[0]
        # a   = coeff[1]
        # E0  = coeff[2]
        # interval = Polyval(zl,z0,zr,dz,b,a,E0,coeff)

        # x = [zl,z0,zr]
        # y = [El,E0,Er]
        # cof = polcof(x,y,2)
        # print('cof =====',cof)
        # print('E0,a,b ',E0,a,b,'\n')
            # END-USE np.polyfit()
    return poly

def Ipoly(z,poly):
    """
    interpolate with polynomial fit
    """
    ix = -1
    for i in range(len(poly)):
        zl = poly[i].zl
        zr = poly[i].zr
        if zl <= z and z <= zr:
            ix = i
            break
    if ix <0:
        raise RuntimeError('arg out of range! {}'.format(z))

            # USE np.polyfit()
            # liefert gleiches Resultat wie Langrange 3 Punkt
    # coeff = poly[ix].coeff
    # res   = np.poly1d(coeff)
    # res   = res(z)
            # END-USE np.poly1d()

    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res

def V0n(poly,n):
        E0 = poly[n].E0                       # [MV/m]
        b  = poly[n].b
        dz = poly[n].dz                       # [cm]
        v0  = E0*(2*dz+2./3.*b*dz**3)*1.e-2   # [MV]
        return v0
def V0(poly,zl2zr):
    v0 = []
    zl,zr = zl2zr
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        v0.append(V0n(poly,i))
    return v0
def Tn(poly,k,n):
        a  = poly[n].a
        b  = poly[n].b
        dz = poly[n].dz
        f1 = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        # DEBUG_MODULE('Tk: (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
        t = f1*f2
        return t
def T(poly,k,zl2zr):
    t = []
    zl,zr = zl2zr
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        # DEBUG_MODULE('Tk: (i,dz,zl,zil,zir,zr)=({:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(i,dz,zl,zil,zir,zr))
        t.append(Tn(poly,k,i))
    return t
def Sn(poly,k,n):
        a  = poly[n].a
        b  = poly[n].b
        dz = poly[n].dz
        s = 2*a*sin(k*dz)/k**2/(2*dz+2./3.*b*dz**3)
        s = s * (1.-k*dz/tan(k*dz))
        return s
def S(poly,k,zl2zr):
    s = []
    zl,zr = zl2zr
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        s.append(Sn(poly,k,i))
    return s
def Tpn(poly,k,n):
        a   = poly[n].a
        b   = poly[n].b
        dz  = poly[n].dz
        tp = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        tp = tp * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        return tp
def Tp(poly,k,zl2zr):
    tp = []
    zl,zr = zl2zr
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        tp.append(Tpn(poly,k,i))
    return tp
def Spn(poly,k,n):
        a   = poly[n].a
        b   = poly[n].b
        dz  = poly[n].dz
        sp = 2*a*sin(k*dz)/k/(2*dz+2./3.*b*dz**2)
        sp = sp * (dz**2-2./k**2+2*dz/k/tan(k*dz))
        return sp
def Sp(poly,k,zl2zr):
    sp = []
    zl,zr = zl2zr
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        sp.append(Spn(poly,k,i))
    return sp

class SFdata(object):
    '''
    Superfish data class  (normiert auf max(E-Feld) = Epeak)
    '''
    def __init__(self,input_file,Epeak=1.):
        print('READING SF-DATA from "{}"'.format(input_file))
        self.input_file = input_file
        self.Epeak      = Epeak
        self.make_Ez_table()
        self.make_Ez_poly()
    def make_Ez_table(self):
        Dpoint = namedtuple('DataPoint',['z','R','Ez'])  # data-point structure
        zp = []
        rp = []
        ep = []
        with open(self.input_file,'r') as f:
            lines = list(f)
            for line in lines[41:-2]:
                # print(line,end='')
                stripped    = line.strip()
                (z,sep,aft) = stripped.partition(' ')
                z =float(z)
                stripped    = aft.strip()
                (R,sep,aft) = stripped.partition(' ')
                R = float(R)
                stripped     = aft.strip()
                (Ez,sep,aft) = stripped.partition(' ')
                Ez = float(Ez)
                zp.append(z)
                rp.append(R)
                ep.append(Ez)
        zprev = [-x for x in reversed(zp[1:])]
        zp = zprev+zp
        rprev = [-x for x in reversed(rp[1:])]
        rp = rprev+rp
        eprev = [x for x in reversed(ep[1:])]
        ep =eprev+ep
        emax = max(ep)
        enorm = self.Epeak/emax
        self.Ez0_tab = []
        for i in range(len(zp)):   # normalize and pack
            ep[i] = ep[i]*enorm
            dpoint = Dpoint(zp[i],rp[i],ep[i])
            self.Ez0_tab.append(dpoint)
    def Ez_table(self):
        return self.Ez0_tab
    def make_Ez_poly(self):
        """
        Calculate polynom coefficients from SF data
        """
        def indexer(nbslices,M):
            """
                nbslices = nboff slices
                N = nboff half-length-slices
                M = nboff SF-points
                n = nboff SF-points/half-length-slice
            """
            N=2*nbslices    # factor 2 more intervals than slices
            if N>M: 
                raise RuntimeError('nboff slices must be <= {}'.format(int((M-1)/2)))
            M = int(M-fmod(M,N))
            n = int(M/N)
            print('{} intervals, {} SF-points, {} SF-points/interval'.format(nbslices,M,2*n))
            for i in range(0,M,2*n):
                # DEBUG_MODULE('(i,i+n,i+2*n)={},{},{}'.format(i,i+n,i+2*n))
                yield((i,i+n,i+2*n))
        
        Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0','coeff'])
        self.poly = []
        anz = 10        # interpolate SF-data with 'anz' polynomials od 2nd order
        for (il,i0,ir) in indexer(anz,len(self.Ez0_tab)):
            # DEBUG_MODULE('(il,i0,ir) ',((il,i0,ir)))
            zl = self.Ez0_tab[il].z
            z0 = self.Ez0_tab[i0].z
            zr = self.Ez0_tab[ir].z
            El = self.Ez0_tab[il].Ez
            E0 = self.Ez0_tab[i0].Ez
            Er = self.Ez0_tab[ir].Ez
            dz = z0-zl
            b  = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
            a  = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
            pval = Polyval(zl,z0,zr,dz,b,a,E0,0.)
            self.poly.append(pval)
    def Ez_poly(self):
        return self.poly

def pre_plt(input_file):
    ax  = plt.subplot(111)
    ax.set_ylabel('Ez0 [MV/m]')
    ax.set_xlabel('z [cm]')
    ax.set_title(input_file)
    return ax

def post_plt(ax):
    plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

def displayLR(table,legend):
    zp   = [+float(x[0]) for x in table]
    Ezp  = [+float(x[2]) for x in table]
    zn   = [-float(x[0]) for x in reversed(table)]
    Ezn  = [+float(x[2]) for x in reversed(table)]
    plt.plot(zn+zp,Ezn+Ezp,label=legend)

def display(table,legend):
    z   = [+float(x[0]) for x in table]
    Ez  = [+float(x[2]) for x in table]
    plt.plot(z,Ez,label=legend)

def test0(input_file):
    Ez0_tab = SFdata(input_file).Ez_table()
    display(Ez0_tab,'SF')
    
def test1():
    '''
    Gauss'che Normalverteilung
    '''
    gap = 4.4
    z = np.arange(0.,gap,gap/500.)
    sigma = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigma,0.)) for x in z]
    displayLR(Ez0_tab,'NG')

            # EXPERIMENTAL
    # EzZ_tab = [(x,0.,Z(x,sigma,0.01 )) for x in z]
    # display(EzZ_tab,'Z')
            # END-EXPERIMENTAL

            # OVERAP?
    # Ez1_tab = [(x,0.,NGauss(x,sigm,gap)) for x in z]
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # cavlen = 2.5
    # print('sigma {}[cm], half-cavity length = {}[cm] = {}sigma'.
    #     format(sigm,cavlen,cavlen/sigm))
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))
            # END-OVERAP?

def test2():
    '''
    Second order polynomial fit with Lagrange 3 point formula to NG formula
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)*1.e-2     # [1/cm]

    anz   = 6            # nboff slices
    gap   = 4.8          # [cm] full gap length
    zl    = -gap/2.      #left  interval boundary
    zr    = gap/2.       #right interval boundary
    sigma = gap/2./1.89  # sigma of NGauss (best fit with SF)
    # sigma = gap/2./2.2   # sigma of NGauss (best fit with SF)
    E0    = 1.           # top of NGauss   (best fit with SF)

    z = np.linspace(zl,zr,2*anz+1)
    Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
    # display(Ez0_tab,'slice')
    poly  = Kpoly(z,sigma,0.,E0*1.)

    zstep = (zr-zl)/500.
    z = np.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez0_tab,'NG-poly')

            # OVERAP?
    # poly1 = Kpoly(z,sigma,gap,E0*1.)  # next cavity
    # Ez1_tab = [(x,0.,Ipoly(x,poly1)) for x in z]
    # display(Ez1_tab,'poly1')
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # display(Ez2_tab,'poly2')
            # END-OVERAP?

    # TTF calculations
    zl2zr = (zl,zr)
    v0  = V0(poly,zl2zr)
    t   = T(poly,k,zl2zr)
    s   = S(poly,k,zl2zr)
    tp  = Tp(poly,k,zl2zr)
    sp  = Sp(poly,k,zl2zr)
    DEBUG_MODULE('V0',v0)
    DEBUG_MODULE('T(k)',t)
    DEBUG_MODULE("T'(k)",tp)
    DEBUG_MODULE('S(k)',s)
    DEBUG_MODULE("S'(k)",sp)

def test3(input_file):
    '''
    Second order polynomial fit with Lagrange 3 point formula to SF data
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)*1.e-2    # [1/cm]
    gap = 4.8
    zl  = -gap/2.
    zr  = +gap/2.
    zl2zr = (zl,zr)

    gap_data = SFdata(input_file,Epeak=1.)
    poly = gap_data.Ez_poly()
    display(gap_data.Ez_table(),'slice')

    zstep = (zr-zl)/500.
    z = np.arange(zl,zr,zstep)
    Ez_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez_tab,'SF-poly')

    # TTF calculations
    v0  = V0(poly,zl2zr)
    t   = T(poly,k,zl2zr)
    s   = S(poly,k,zl2zr)
    tp  = Tp(poly,k,zl2zr)
    sp  = Sp(poly,k,zl2zr)
    DEBUG_MODULE('V0',v0)
    DEBUG_MODULE('T(k)',t)
    DEBUG_MODULE("T'(k)",tp)
    DEBUG_MODULE('S(k)',s)
    DEBUG_MODULE("S'(k)",sp)

if __name__ == '__main__':
    input_file='SF_WDK2g44.TBL'
    ax = pre_plt(input_file)
    test0(input_file)               # SF
    test1()                         # NG
    test2()                         # poly-fit to NG
    test3(input_file)               # poly-fit to SF
    post_plt(ax)
