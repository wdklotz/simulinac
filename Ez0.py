#!/Users/klotz/anaconda3/bin/python3.6
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
import sys
import matplotlib.pyplot as plt
import numpy as NP
from math import sin,tan,pi,exp,fmod,cos
from collections import namedtuple
import pprint, inspect
import unittest

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')
DEBUG_MODULE = DEBUG_OFF
DEBUG_TEST2  = DEBUG_OFF
DEBUG_TEST3  = DEBUG_OFF

from setutil import PARAMS,Proton
# Polyval: polynomial approximation for E(z,r=0), z in interval [zl,zr]: see (4.4.1) A.Shishlo/J.Holmes
Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0','coeff'])
# Dpoint: Table data point -  _Ez0_tab is list(Dpoint)
Dpoint  = namedtuple('Dpoint',['z','R','Ez'])

def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res
def Kpoly(z,sigma,mu,E):
    """Calculate polynom coefficients from NG formula"""
    poly = []
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
        a = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
        pval = Polyval(zl,z0,zr,dz,b,a,E0,0.)
        poly.append(pval)

            # USING NP.polyfit() liefert gleiches Resultat wie Langrange 3 Punkt
        # x   = NP.array((zl,z0,zr))
        # y   = NP.array((El,E0,Er))
        # coeff = NP.polyfit(x,y,2)
        # b   = coeff[0]
        # a   = coeff[1]
        # E0  = coeff[2]
        # interval = Polyval(zl,z0,zr,dz,b,a,E0,coeff)

        # x = [zl,z0,zr]
        # y = [El,E0,Er]
        # cof = polcof(x,y,2)
        # print('cof =====',cof)
        # print('E0,a,b ',E0,a,b,'\n')
            # END-USING NP.polyfit()
    return poly
def Ipoly(z,poly):
    """Interpolation using the polynomial fit"""
    ix = -1
    for i in range(len(poly)):
        zl = poly[i].zl
        zr = poly[i].zr
        if zl <= z and z <= zr:
            ix = i
            break
    if ix <0:
        raise RuntimeError('Ipoly(): arg out of range! {}'.format(z))

    #       USE NP.polyfit()
    # liefert gleiches Resultat wie Langrange 3 Punkt
    # coeff = poly[ix].coeff
    # res   = NP.poly1d(coeff)
    # res   = res(z)
    #       END-USE NP.poly1d()

    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res
def V0n(poly,n):
    """ Formel (4.4.3) A.Shishlo/J.Holmes """
    E0 = poly[n].E0                       # [MV/m]
    b  = poly[n].b
    dz = poly[n].dz                       # [cm]
    v0  = E0*(2*dz+2./3.*b*dz**3)*1.e-2   # [MV]
    return v0
def V0(poly,zintval):
    v0 = []
    zl,zr = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        if zil < zl or zir > zr: continue
        v0.append(V0n(poly,i))
    return v0
def Tn(poly,k,n):
    """ Formel (4.4.6) A.Shishlo/J.Holmes """
    a  = poly[n].a
    b  = poly[n].b
    dz = poly[n].dz
    f1 = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
    f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
    DEBUG_MODULE('Tn(): (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
    t = f1*f2
    return t
def T(poly,k,zintval):
    t = []
    zl,zr = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz = poly[i].dz
        if zil < zl or zir > zr: continue
        DEBUG_MODULE('T(): (i,dz,zl,zil,zir,zr)=({:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(i,dz,zl,zil,zir,zr))
        t.append(Tn(poly,k,i))
    return t
def Sn(poly,k,n):
    """ Formel (4.4.7) A.Shihlo """
    a  = poly[n].a
    b  = poly[n].b
    dz = poly[n].dz
    s = 2*a*sin(k*dz)/k**2/(2*dz+2./3.*b*dz**3)
    s = s * (1.-k*dz/tan(k*dz))
    return s
def S(poly,k,zintval):
    s = []
    zl,zr = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        if zil < zl or zir > zr: continue
        s.append(Sn(poly,k,i))
    return s
def Tpn(poly,k,n):
    """ Formel (4.4.8) A.Shishlo/J.Holmes """
    b   = poly[n].b
    dz  = poly[n].dz
    tp = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
    tp = tp * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
    return tp
def Tp(poly,k,zintval):
    tp = []
    zl,zr = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        if zil < zl or zir > zr: continue
        tp.append(Tpn(poly,k,i))
    return tp
def Spn(poly,k,n):
    """ Formel (4.4.9) A.Shishlo/J.Holmes """
    a   = poly[n].a
    b   = poly[n].b
    dz  = poly[n].dz
    sp = 2*a*sin(k*dz)/k/(2*dz+2./3.*b*dz**2)
    sp = sp * (dz**2-2./k**2+2*dz/k/tan(k*dz))
    return sp
def Sp(poly,k,zintval):
    sp = []
    zl,zr = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        if zil < zl or zir > zr: continue
        sp.append(Spn(poly,k,i))
    return sp
class SFdata(object):
    ''' Cavity E(z,r=0) field profile: Superfish data  (can be normiert auf EzPeak & gap)
    IN: EzPeak peak field [MV/m]
    gap: full gap [m]
    '''
    def __init__(self,input_file,EzPeak=0.,gap=0.):
        print('READING SF-DATA from "{}"'.format(input_file))
        self.input_file = input_file
        self._EzPeak    = EzPeak
        self._gap       = gap           # full gap !!
        self._EzAvg     = None          # for later
        self._Ez0_tab_raw,EzAvg = self.make_Ez_table()  # rtaw data from SF will never be scaled!
        self._Ez0_tab_scaled,self._EzAvg = self.scale_Ez_table(EzAvg,EzPeak=self._EzPeak,gap=self._gap) # scale EzAvg and profile
        self._poly = self.make_polyfit(self._Ez0_tab_scaled)

    @property
    def gap(self):
        return self._gap
    @property
    def EzAvg(self):
        return self._EzAvg
    @property
    def EzPeak(self):
        return self._EzPeak
    @property
    def EzPoly(self):
        return self._poly
    def make_Ez_table(self):
        """ read raw data and scale to self.EzPeak and self.gap """
        zp = []; rp = []; ep = []
        leading  =  41      # nbof leading lines to skip
        trailing  = 2       # nbof trailing lines to skip
        raw_tab = []
        with open(self.input_file,'r') as f:
            lines = list(f)
            lines = lines[leading:-trailing]           # remove leading and trailing lines

            for line in lines:
                stripped    = line.strip()
                DEBUG_OFF(stripped)
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
        f.close()
        zprev = [-x for x in reversed(zp[1:])]
        zp = zprev+zp
        rprev = [x for x in reversed(rp[1:])]
        rp = rprev+rp
        eprev = [x for x in reversed(ep[1:])]
        ep = eprev+ep
        N = len(zp)
        EzAvg = 0
        for i in range(N):  
            raw_tab.append(Dpoint(zp[i],rp[i],ep[i]))
            EzAvg += ep[i]
        EzAvg = EzAvg/N
        return raw_tab,EzAvg
    def scale_Ez_table(self,EzAvg,EzPeak=0.,gap=0.):   # NOTE: full gap!
        # scales Ez-axis and z-axis of raw SF-data-table
        EzMax = max([x.Ez for x in self._Ez0_tab_raw])
        self._EzPeak = EzMax if EzPeak == 0. else EzPeak
        zmax  = max([x.z  for x in self._Ez0_tab_raw])
        self._gap = 2.*zmax if gap == 0. else gap

        tab = []
        half_gap = self._gap*0.5
        for i in range(len(self._Ez0_tab_raw)):  
            z = self._Ez0_tab_raw[i].z*(half_gap/zmax)       # scale z-axis
            e = self._Ez0_tab_raw[i].Ez*(self._EzPeak/EzMax)  # sclae Ez-axis
            r = self._Ez0_tab_raw[i].R*(half_gap/zmax)       # scale R same as Z
            tab.append(Dpoint(z,r,e))

        EzAvg = EzAvg*(self._EzPeak/EzMax)
        return tab,EzAvg
    def make_polyfit(self,Ez_table): 
        # Polynomial fits to raw data according to Shislo&Holmes
        # In here M adjacent raw intervals are taken as one single interval and fitted with a polynomial
        sf_tab = Ez_table        # SF-table scales applied
        N      = len(sf_tab)     # stuetzpunkte in raw table
        M      = 8               # raw intervalls/poly interval (must be even number [2,4,6,8,....])
        polies = []              # PolyFit: list(Polyvals)

        DEBUG_HERE = DEBUG_OFF

        DEBUG_HERE('Ez0_poly::SFdata::make_polyfit: raw function values: {} in {}'.format(N,range(N-1)))
        DEBUG_HERE('Ez0_poly::SFdata::make_polyfit: first is sf_tab[{:3}]..{}'.format(0,sf_tab[0]))
        DEBUG_HERE('Ez0_poly::SFdata::make_polyfit: last is  sf_tab[{:3}]..{}'.format(N-1,sf_tab[N-1]))

        i=0
        while(True):
            il = i
            i0 = i+int(M/2)
            ir = i+M
            if i>N-1 or i0>N-1 or ir>N-1: break
            zl = sf_tab[il].z
            z0 = sf_tab[i0].z
            zr = sf_tab[ir].z
            El = sf_tab[il].Ez
            E0 = sf_tab[i0].Ez
            Er = sf_tab[ir].Ez
            i = ir   # next interval
            dz = z0-zl
            b  = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
            a  = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
            pval = Polyval(zl,z0,zr,dz,b,a,E0,0.)
            polies.append(pval)
            DEBUG_HERE('Ez0_poly::SFdata::make_polyfit: (il,i0,ir) ({:3},{:3},{:3}),  (zl,z0,zr,E0) ({:6.3f},{:6.3f},{:6.3f},{:6.3f})'.format(il,i0,ir,zl,z0,zr,E0))
        print('Ez0_poly::SFdata::make_polyfit: {} poly intervals'.format(len(polies)))
        # self._poly = polies
        return polies
    def Ez0t(self, z, t, omega, phis):
        """E(z,0,t): time dependent field value at location z"""
        res = Ipoly(z,self.EzPoly) * cos(omega*t+phis)
        return res
    def dEz0tdt(self, z, t, omega, phis):
        """dE(z,0,t)/dt: time derivative of field value at location z"""
        res = - omega * Ipoly(z,self.EzPoly) * sin(omega*t+phis)
        return res
class TestEz0Methods(unittest.TestCase):
    @classmethod
    def displayLR(cls,table,legend):
        # display left (L) and right (R) table data
        zp   = [+float(x[0]) for x in table]
        Ezp  = [+float(x[2]) for x in table]
        zn   = [-float(x[0]) for x in reversed(table)]
        Ezn  = [+float(x[2]) for x in reversed(table)]
        plt.plot(zn+zp,Ezn+Ezp,label=legend)
    @classmethod
    def display(cls,table,legend):
        z   = [+float(x[0]) for x in table]
        Ez  = [+float(x[2]) for x in table]
        plt.plot(z,Ez,label=legend)
    def test_0(self):
        print("\b----------------------------------------test_0")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        sfdata  = SFdata(input_file,EzPeak=1.5,gap=5.)   #NOTE: full gap
        # sfdata  = SFdata(input_file)
        Ez0_tab        = sfdata._Ez0_tab_raw
        Ez1_tab        = sfdata._Ez0_tab_scaled
        Ez2_tab, dummy = sfdata.scale_Ez_table(sfdata._EzAvg,EzPeak=0.8,gap=6.)

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'SF-raw')
        self.display(Ez1_tab,'SF-scaled.1')
        self.display(Ez2_tab,'SF-scaled.2')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(input_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test_1(self):
        print("\b----------------------------------------test_1")
        gap = 4.4
        z = NP.arange(0.,gap,gap/500.)
        sigma = 1.14
        Ez0_tab = [(x,0.,NGauss(x,sigma,0.)) for x in z]

        ax  = plt.subplot(111)
        self.displayLR(Ez0_tab,'NG')    
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test_2(self):
        print("\b----------------------------------------test_2")
        particle = Proton(100.)
        beta     = particle.beta
        c        = PARAMS['clight']
        freq     = 800.e6
        k        = 2*pi*freq/(c*beta)*1.e-2     # [1/cm]

        anz   = 6            # nboff slices
        gap   = 4.8          # [cm] full gap length
        zl    = -gap/2.      #left  interval boundary
        zr    = gap/2.       #right interval boundary
        sigma = gap/2./1.89  # sigma of NGauss (best fit with SF)
        # sigma = gap/2./2.2   # sigma of NGauss (best fit with SF)
        E0    = 1.           # top of NGauss   (best fit with SF)

        z = NP.linspace(zl,zr,2*anz+1)
        Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
        # display(Ez0_tab,'slice')
        poly  = Kpoly(z,sigma,0.,E0*1.)

        zstep = (zr-zl)/500.
        z = NP.arange(zl,zr,zstep)
        Ez0_tab = [(x,0.,Ipoly(x, poly)) for x in z]

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'NG-poly')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test_3(self):
        print("\b----------------------------------------test_3")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        EzPeak = 1.2
        gap    = 0.
        sf_data   = SFdata(input_file,EzPeak=EzPeak,gap=gap)
        poly_data = sf_data._poly
        particle  = Proton(tkin=100.)
        beta      = particle.beta
        c         = PARAMS['clight']
        freq      = 800.e6
        k         = 2*pi*freq/(c*beta)*1.e-2    # [1/cm]
        zl        = -sf_data.gap/2.
        zr        = -zl
        zintval   = (zl,zr)
        zstep     = (zr-zl)/500.

        z = NP.arange(zl,zr,zstep)
        ipoly_werte = [(x,0.,Ipoly(x, poly_data)) for x in z]

        ax  = plt.subplot(111)
        self.display(sf_data._Ez0_tab_scaled,'SFdata')
        self.display(ipoly_werte,'Polydata')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(input_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
        # TTF calculations
        v0  = V0(poly_data,zintval)
        t   = T(poly_data,k,zintval)
        s   = S(poly_data,k,zintval)
        tp  = Tp(poly_data,k,zintval)
        sp  = Sp(poly_data,k,zintval)
        DEBUG_TEST3('V0 {}'.format(v0))
        DEBUG_TEST3('T(k) {}'.format(t))
        DEBUG_TEST3("T'(k) {}".format(tp))
        DEBUG_TEST3('S(k) {}'.format(s))
        DEBUG_TEST3("S'(k) {}".format(sp))
    def test_4(self):
        print("\b----------------------------------------test_4")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        EzPeak = 1.2
        gap    = 0.
        sfdata = SFdata(input_file, EzPeak=EzPeak,gap=gap)
        print("peak:{} -- average:{} -- average/peak {}".format(sfdata.EzPeak,sfdata.EzAvg,sfdata.EzAvg/sfdata.EzPeak))

if __name__ == '__main__':
     unittest.main()   
