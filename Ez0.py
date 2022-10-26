#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v10.22.7'
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
import numpy as NP
import unittest
from operator import eq
from math import sin,tan,pi,exp,fmod,cos
from collections import namedtuple

from setutil import PARAMS,Proton,DEBUG_ON,DEBUG_OFF
# Polyval: polynomial approximation for E(z,r=0), z in interval [zl,zr]: see (4.4.1) A.Shishlo/J.Holmes
Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0','coeff'])
# Dpoint: Table data point -  Ez0_tab is list(Dpoint)
Dpoint  = namedtuple('Dpoint',['z','R','Ez'])

def NGauss(x,sig,mu):
     # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res
def Kpoly(z,sigma,mu,E):
    #  Polynom coefficients from Gauss Normalverteilung
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
    if ix < 0:
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
    DEBUG_OFF('Tn(): (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
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
        DEBUG_OFF('T(): (i,dz,zl,zil,zir,zr)=({:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(i,dz,zl,zil,zir,zr))
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
    """ 
    Cavity E(z,r=0) field profile from Superfish.
    Raw data can be scaled to fit EzPeak & gap
    SFdata.field_data(input_file,EzPeak,gap) will return a singleton SFdata object 
        bound to input_file and scaled to EzPeak and gap.
    IN: EzPeak peak field [MV/m], gap: full gap [m]
    OUT: SFdata object
    """
    # class variable
    instances = {}    
    @classmethod
    def field_data(cls,input_file,EzPeak=0.,gap=0.):
        def scale(instance,EzPeak,gap):
            instance.EzPeak = EzPeak
            instance.gap    = gap
            instance.Ez0_tab, instance.EzAvg = instance.scale_Ez_table(EzPeak,gap)    # scale EzAvg and profile
            instance.poly = instance.make_polyfit(instance.Ez0_tab)
            DEBUG_OFF('make_polyfit: {} poly intervals'.format(len(instance.poly)))
        instance = cls.instances.get(input_file)
        if instance == None:
            # new instance
            instance = SFdata(input_file)
            if EzPeak != 0. and gap != 0.:
                scale(instance,EzPeak,gap)
            elif EzPeak != 0. and gap == 0.:
                scale(instance,EzPeak,instance.gap_raw)
            elif EzPeak == 0. and gap != 0.:
                scale(instance,instance.EzPeak_raw,gap)
            elif EzPeak == 0. and gap == 0.:
                pass
            SFdata.instances[input_file] = instance
        else:
            # return [scaled] instance 
            if EzPeak != 0. and gap != 0. and (EzPeak != instance.EzPeak or gap != instance.gap):
                scale(instance,EzPeak,gap)
            elif EzPeak != 0. and gap == 0.:
                scale(instance,EzPeak,instance.gap)
            elif EzPeak == 0. and gap != 0.:
                scale(instance,instance.EzPeak,gap)
            elif EzPeak == 0. and gap == 0. or (EzPeak == instance.EzPeak and gap == instance.gap):
                pass
        return instance
    @property
    def EzPoly(self):
        return self.poly

    def __init__(self,input_file):
        print('READING SF-DATA from "{}"'.format(input_file))
        self.input_file = input_file
        self.Ez0_tab_raw, self.EzAvg_raw = self.make_Ez_table(input_file)  # raw data from SF will never be modified!
        self.EzPeak_raw = max([x.Ez for x in self.Ez0_tab_raw])
        self.gap_raw    = 2.*max([x.z  for x in self.Ez0_tab_raw])

        self.Ez0_tab    = self.Ez0_tab_raw
        self.EzAvg      = self.EzAvg_raw
        self.EzPeak     = self.EzPeak_raw
        self.gap        = self.gap_raw
        self.poly       = self.make_polyfit(self.Ez0_tab)

    def make_Ez_table(self,input_file):
        """ read raw data and scale to self.EzPeak and self.gap """
        zp = []; rp = []; ep = []
        leading  =  41      # nbof leading lines to skip
        trailing  = 2       # nbof trailing lines to skip
        raw_tab = []
        with open(input_file,'r') as f:
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
    def scale_Ez_table(self,EzPeak,gap):   # NOTE: full gap!
        Ez0_tab = []
        # half_gap = self.gap_raw*0.5
        EzMax = self.EzPeak_raw
        zmax  = self.gap_raw
        for i in range(len(self.Ez0_tab_raw)):
            z = self.Ez0_tab_raw[i].z*(gap/zmax)       # scale z-axis
            e = self.Ez0_tab_raw[i].Ez*(EzPeak/EzMax)  # sclae Ez-axis
            r = self.Ez0_tab_raw[i].R*(gap/zmax)       # scale R same as Z
            Ez0_tab.append(Dpoint(z,r,e))

        # self.poly = self.make_polyfit(Ez0_tab)
        EzAvg     = self.EzAvg_raw*(EzPeak/EzMax)
        return Ez0_tab,EzAvg
    def make_polyfit(self,Ez_table):
        # Polynomial fits to raw data according to Shislo&Holmes
        # In here M adjacent raw intervals are taken as one single interval and fitted with a polynomial
        sf_tab = Ez_table        # SF-table scales applied
        N      = len(sf_tab)     # stuetzpunkte in raw table
        M      = 8               # raw intervalls/poly interval (must be even number [2,4,6,8,....])
        polies = []              # PolyFit: list(Polyvals)

        DEBUG_OFF('make_polyfit: raw function values: {} in {}'.format(N,range(N-1)))
        DEBUG_OFF('make_polyfit: first is sf_tab[{:3}]..{}'.format(0,sf_tab[0]))
        DEBUG_OFF('make_polyfit: last is  sf_tab[{:3}]..{}'.format(N-1,sf_tab[N-1]))
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
            DEBUG_OFF('Ez0_poly::SFdata::make_polyfit: (il,i0,ir) ({:3},{:3},{:3}),  (zl,z0,zr,E0) ({:6.3f},{:6.3f},{:6.3f},{:6.3f})'.format(il,i0,ir,zl,z0,zr,E0))
        DEBUG_OFF('make_polyfit: {} poly intervals'.format(len(polies)))
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
    def test0(self):
        print("\b---------------------------------test0 (singleton)")
        input_file='SF/SF_WDK2g44.TBL'
        (EzPeak,gap) = (4.0,4.0)
        sfdata1 = SFdata.field_data(input_file,*(EzPeak,gap))
        print('   ',(EzPeak,gap),sfdata1.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        sfdata2 = SFdata.field_data(input_file,*(EzPeak,gap))
        print('   ',(EzPeak,gap),sfdata2.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,gap) = (5.0,2.0)
        sfdata3 = SFdata.field_data(input_file,*(EzPeak,gap))
        print('   ',(EzPeak,gap),sfdata3.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,gap) = (0.,0.)
        sfdata4 = SFdata.field_data(input_file)
        print('   ',(EzPeak,gap),sfdata4.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        input_file='SF/SF_WDK2g22.TBL'
        sfdata5 = SFdata.field_data(input_file)
        print('   ',(EzPeak,gap),sfdata5.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,gap) = (5.0,2.0)
        sfdata6 = SFdata.field_data(input_file,*(EzPeak,gap))
        print('   ',(EzPeak,gap),sfdata6.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        input_file='SF/SF_WDK2g44.TBL'
        (EzPeak,gap) = (0.,0.)
        sfdata7 = SFdata.field_data(input_file)
        print('   ',(EzPeak,gap),sfdata7.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])
    def test1(self):
        print("\b----------------------------------------test1")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        sfdata  = SFdata.field_data(input_file,EzPeak=1.5,gap=5.)   #NOTE: full gap
        # sfdata  = SFdata(input_file)
        Ez0_tab        = sfdata.Ez0_tab_raw
        Ez1_tab        = sfdata.Ez0_tab
        Ez2_tab, dummy = sfdata.scale_Ez_table(EzPeak=2.5,gap=2.)

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'SF-raw')
        self.display(Ez1_tab,'SF-scaled.1')
        self.display(Ez2_tab,'SF-scaled.2')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(input_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test2(self):
        print("\b----------------------------------------test2")
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
    def test3(self):
        print("\b----------------------------------------test3")
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
    def test4(self):
        print("\b----------------------------------------test4")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        EzPeak = 1.2
        gap    = 0.
        sf_data   = SFdata.field_data(input_file,EzPeak=EzPeak,gap=0.)
        poly_data = sf_data.poly
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
        self.display(sf_data.Ez0_tab,'SFdata')
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
        DEBUG_OFF('V0 {}'.format(v0))
        DEBUG_OFF('T(k) {}'.format(t))
        DEBUG_OFF("T'(k) {}".format(tp))
        DEBUG_OFF('S(k) {}'.format(s))
        DEBUG_OFF("S'(k) {}".format(sp))
    def test5(self):
        print("\b----------------------------------------test5")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        EzPeak = 1.2
        gap    = 0.
        sfdata = SFdata.field_data(input_file,EzPeak=EzPeak,gap=gap)
        print("peak:{} -- average:{} -- average/peak {}".format(sfdata.EzPeak,sfdata.EzAvg,sfdata.EzAvg/sfdata.EzPeak))

if __name__ == '__main__':
    # unittest.main()
    tests = TestEz0Methods()
    tests.test0()    
    tests.test1()
    tests.test2()
    tests.test3()
    tests.test4()
    tests.test5()
    