#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='11.0.2.3'
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
Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0'])
# Dpoint: Table data point - Ez0_tab is list(Dpoint)
Dpoint  = namedtuple('Dpoint',['z','R','Ez'])

def NormalGauss(x,sig,mu):
     # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res
def GaussPoly(z,sigma,mu,E):
    #  Polynom coefficients from Gauss Normalverteilung
    poly = []
    anz = len(z)
    for i in range(0,anz-2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  = z0-zl
        El  = E*NormalGauss(zl,sigma,mu)
        E0  = E*NormalGauss(z0,sigma,mu)
        Er  = E*NormalGauss(zr,sigma,mu)
        b = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation
        a = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
        pval = Polyval(zl,z0,zr,dz,b,a,E0)
        poly.append(pval)
    return poly
def EzPoly(z,poly):
    """ Interpolation using the polynomial fit """
    ix = -1
    for i in range(len(poly)):
        zl = poly[i].zl
        zr = poly[i].zr
        if zl <= z and z <= zr:
            ix = i
            break
    if ix < 0:
        raise RuntimeError('EzPoly(): arg out of range! {}'.format(z))

    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res
def zEPoly(EzAvg,polyValues):
    """ Interpolation of Umkehr-function z(E) using a polynomial fit to get zeff=z(EzAvg) for argument EzAvg.
    IN:
        EzAvg - energy average calculated from SF distribution
        polyValues - the multi-interval polynomial fits to the the scaled SF distribution
    OUT:
        zeff = z(EzAvg) which is effictive 1/2 gap size
    """
    ix=-1              # search from right to left (positive z values) for intervals containingin vicinity of EzAvg
    bwhile = True
    while(bwhile):
        polyval = polyValues[ix]
        e0 = polyval.E0
        # zl = polyval.zl
        # zr = polyval.zr
        # dz = polyval.dz
        # DEBUG_ON(f"[ix={ix},zl={zl},zr={zr}],dz={dz},E0={e0}") 
        if e0 < EzAvg:
            ix = ix-1
        else:
            bwhile = False
    # prepare for numpy's polyfit over two adjacent intervals, one of them containing EzAvg
    z2 = polyValues[ix+1].zr
    z1 = polyValues[ix+1].zl
    z0 = polyValues[ix].zl
    e2 = EzPoly(z2,polyValues)
    e1 = EzPoly(z1,polyValues)
    e0 = EzPoly(z0,polyValues)
    DEBUG_OFF(f"z0 {z0},z1 {z1},z2 {z2},e0 {e0},e1 {e1},e2 {e2},")
    #   use NP.polyfit(...,deg=2)
    pfit = NP.polyfit([e0,e1,e2],[z0,z1,z2],deg=2)
    zeff = pfit[0]*EzAvg**2 + pfit[1]*EzAvg + pfit[2]
    DEBUG_ON(f"input EzAvg {EzAvg}, fitted EzAvg(zeff) {EzPoly(zeff,polyValues)}, zeff {zeff}")
    return zeff
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
    Cavity E(z,r=0) a.k.a. Ez0 field profile from Superfish xxx.TBL file.
    Raw data can be scaled on Ordinate and Abscisse to fit EzPeak & definition interval.
    SFdata.field_data(input_file,EzPeak,IntgInterval) will return a singleton SFdata object bound to input_file.
    IN: EzPeak peak field [MV/m]
        IntgInterval: full interval [m] on which Ez0 is defined.
    OUT: SFdata object
    """
    # class variables
    instances = {}    # dict of sigleton objects
    @classmethod
    def field_data(cls,input_file,EzPeak=0.,IntgInterval=0.):
        def scale(sfdata,EzPeak,IntgInterval):
            sfdata.EzPeak          = EzPeak
            sfdata.IntgInterval    = IntgInterval
            (sfdata.Ez0_tab, sfdata.EzAvg) = sfdata.scale_Ez_table(EzPeak,IntgInterval)    #NOTE: IntgInterval is length [-z,+z]
            # fit polynomials after scaling
            sfdata._Ezpoly = sfdata.make_Ezpoly(sfdata.Ez0_tab)
            DEBUG_OFF('make_Ezpoly: {} poly intervals'.format(len(sfdata.polyValues)))

        """ field_data body -------- field_data body -------- field_data body -------- field_data body -------- field_data body """
        """ field_data body -------- field_data body -------- field_data body -------- field_data body -------- field_data body """
        # use a key made from <SFdata file name>%<Ezpeak value>%<IntgInterval value> to retrieve the corresponding field table instance
        instance_key = '{}%{}%{}'.format(input_file,EzPeak,IntgInterval)
        # instance with instance_key exists? get it; else None
        instance = cls.instances.get(instance_key)
        if instance == None:
            # new scaled instance
            instance = SFdata(input_file)
            if EzPeak != 0. and IntgInterval != 0.:
                scale(instance,EzPeak,IntgInterval)
            elif EzPeak != 0. and IntgInterval == 0.:
                scale(instance,EzPeak,instance.IntgInterval_raw)
            elif EzPeak == 0. and IntgInterval != 0.:
                scale(instance,instance.EzPeak_raw,IntgInterval)
            elif EzPeak == 0. and IntgInterval == 0.:
                scale(instance,instance.EzPeak_raw,instance.IntgInterval_raw)
            SFdata.instances[instance_key] = instance
            DEBUG_OFF(SFdata.instances)
        return instance
    @property
    def polyValues(self):
        return self._Ezpoly
    def Ez0t(self, z, t, omega, phis):
        """E(z,0,t): time dependent field value at location z"""
        res = EzPoly(z,self.polyValues) * cos(omega*t+phis)
        return res
    def dEz0tdt(self, z, t, omega, phis):
        """dE(z,0,t)/dt: time derivative of field value at location z"""
        res = - omega * EzPoly(z,self.polyValues) * sin(omega*t+phis)
        return res

    def __init__(self,input_file):
        print('READING SF-DATA from "{}"'.format(input_file))
        self.input_file = input_file
        self.make_Ez_table(input_file)  # Ez0_tab_raw is full table on symetric interval [-z,+z]
        # actual [scaled] data will go here
        self.Ez0_tab    = self.Ez0_tab_raw
        self.EzAvg      = self.EzAvg_raw
        self._Ezpoly      = None  # will be set when scaled
    def make_Ez_table(self,input_file):
        """ read raw data and return raw table and raw table's average """
        zp = []; rp = []; ep = []
        leading  =  41      # nbof leading lines to skip
        trailing  = 2       # nbof trailing lines to skip
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

        # construct table for negative z values (mirror on vertical axis)
        zprev = [-x for x in reversed(zp[1:])]
        zp = zprev+zp
        rprev = [x for x in reversed(rp[1:])]
        rp = rprev+rp
        eprev = [x for x in reversed(ep[1:])]
        ep = eprev+ep
        N = len(zp)
        EzAvg = 0
        raw_tab = []
        for i in range(N):
            raw_tab.append(Dpoint(zp[i],rp[i],ep[i]))
            EzAvg += ep[i]
        EzAvg = EzAvg/N
        # raw data from SF will never be modified!
        self.Ez0_tab_raw       = raw_tab
        self.EzAvg_raw         = EzAvg
        self.EzPeak_raw        = max([x.Ez for x in self.Ez0_tab_raw])
        self.IntgInterval_raw  = 2.*self.Ez0_tab_raw[-1].z   #NOTE: IntgInterval_raw is length [-z,+z]
        return
    def scale_Ez_table(self,EzPeak,IntgInterval):   #NOTE: IntgInterval is length [-z,+z]
        Ez0_tab = []
        EzMax = self.EzPeak_raw
        zmax  = self.IntgInterval_raw
        for i in range(len(self.Ez0_tab_raw)):
            z = self.Ez0_tab_raw[i].z*(IntgInterval/zmax)       # scale z-axis
            e = self.Ez0_tab_raw[i].Ez*(EzPeak/EzMax)           # sclae Ez-axis
            r = self.Ez0_tab_raw[i].R*(IntgInterval/zmax)       # scale R same as Z
            Ez0_tab.append(Dpoint(z,r,e))
        EzAvg     = self.EzAvg_raw*(EzPeak/EzMax)
        return Ez0_tab,EzAvg
    def make_Ezpoly(self,Ez_table):
        # Polynomial fits to raw data according to Shislo&Holmes
        # In here M adjacent raw intervals are taken as one single interval and fitted with a polynomial
        sf_tab = Ez_table        # SF-table scales applied
        N      = len(sf_tab)     # stuetzpunkte in raw table
        M      = 8               # raw intervalls/poly interval (must be even number [2,4,6,8,....])
        polies = []              # PolyFit: list(Polyvals)

        DEBUG_OFF('make_Ezpoly: raw function values: {} in {}'.format(N,range(N-1)))
        DEBUG_OFF('make_Ezpoly: first is sf_tab[{:3}]..{}'.format(0,sf_tab[0]))
        DEBUG_OFF('make_Ezpoly: last is  sf_tab[{:3}]..{}'.format(N-1,sf_tab[N-1]))
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
            pval = Polyval(zl,z0,zr,dz,b,a,E0)
            polies.append(pval)
            DEBUG_OFF('Ez0_poly::SFdata::make_Ezpoly: (il,i0,ir) ({:3},{:3},{:3}),  (zl,z0,zr,E0) ({:6.3f},{:6.3f},{:6.3f},{:6.3f})'.format(il,i0,ir,zl,z0,zr,E0))
        DEBUG_OFF('make_Ezpoly: {} poly intervals'.format(len(polies)))
        return polies

class TestEz0Methods(unittest.TestCase):
    @classmethod
    def display_with_mirror(cls,table,legend):
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
        (EzPeak,IntgInterval) = (4.0,4.0)
        sfdata1 = SFdata.field_data(input_file,*(EzPeak,IntgInterval))
        print('   ',(EzPeak,IntgInterval),sfdata1.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        sfdata2 = SFdata.field_data(input_file,*(EzPeak,IntgInterval))
        print('   ',(EzPeak,IntgInterval),sfdata2.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,IntgInterval) = (5.0,2.0)
        sfdata3 = SFdata.field_data(input_file,*(EzPeak,IntgInterval))
        print('   ',(EzPeak,IntgInterval),sfdata3.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,IntgInterval) = (0.,0.)
        sfdata4 = SFdata.field_data(input_file)
        print('   ',(EzPeak,IntgInterval),sfdata4.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        input_file='SF/SF_WDK2g22.TBL'
        sfdata5 = SFdata.field_data(input_file)
        print('   ',(EzPeak,IntgInterval),sfdata5.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,IntgInterval) = (5.0,2.0)
        sfdata6 = SFdata.field_data(input_file,*(EzPeak,IntgInterval))
        print('   ',(EzPeak,IntgInterval),sfdata6.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        input_file='SF/SF_WDK2g44.TBL'
        (EzPeak,IntgInterval) = (0.,0.)
        sfdata7 = SFdata.field_data(input_file)
        print('   ',(EzPeak,IntgInterval),sfdata7.input_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])
    def test1(self):
        print("\b----------------------------------------test1")
        # input_file='SF/PILL-2CM.TBL'
        # input_file='SF/SF_WDK2g44.TBL'
        input_file='SF/CAV-FLAT-R135-L31.TBL'
        sfdata  = SFdata.field_data(input_file,EzPeak=1.5,IntgInterval=5.)   #NOTE: IntgInterval is length [-z,+z]
        Ez0_tab        = sfdata.Ez0_tab_raw
        Ez1_tab        = sfdata.Ez0_tab
        Ez2_tab, dummy = sfdata.scale_Ez_table(EzPeak=2.5,IntgInterval=2.)   #NOTE: IntgInterval is length [-z,+z]

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
        IntgInterval = 4.4
        z = NP.arange(0.,IntgInterval,IntgInterval/500.)
        sigma = 1.14
        Ez0_tab = [(x,0.,NormalGauss(x,sigma,0.)) for x in z]

        ax  = plt.subplot(111)
        self.display_with_mirror(Ez0_tab,'NG')
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
        IntgInterval = 4.8          # [cm] full interval length
        zl    = -IntgInterval/2.      #left  interval boundary
        zr    = IntgInterval/2.       #right interval boundary
        sigma = IntgInterval/2./1.89  # sigma of NormalGauss (best fit with SF)
        # sigma = IntgInterval/2./2.2   # sigma of NormalGauss (best fit with SF)
        E0    = 1.           # top of NormalGauss   (best fit with SF)

        z = NP.linspace(zl,zr,2*anz+1)
        Ez0_tab = [(x,0.,E0*NormalGauss(x,sigma,0.)) for x in z]
        # display(Ez0_tab,'slice')
        poly  = GaussPoly(z,sigma,0.,E0*1.)

        zstep = (zr-zl)/500.
        z = NP.arange(zl,zr,zstep)
        Ez0_tab = [(x,0.,EzPoly(x, poly)) for x in z]

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'NG-poly')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test4(self):
        print("\b----------------------------------------test4")
        # input_file='SF/PILL-2CM.TBL'
        # input_file='SF/SF_WDK2g44.TBL'
        input_file='SF/CAV-FLAT-R135-L31.TBL'
        EzPeak = 1.2
        IntgInterval    = 0.
        sf_data   = SFdata.field_data(input_file,EzPeak=EzPeak,IntgInterval=0.)
        polyValues= sf_data.polyValues
        particle  = Proton(tkin=100.)
        beta      = particle.beta
        c         = PARAMS['clight']
        freq      = 800.e6
        k         = 2*pi*freq/(c*beta)*1.e-2    # [1/cm]
        zl        = -sf_data.IntgInterval/2.
        zr        = -zl
        zintval   = (zl,zr)
        zstep     = (zr-zl)/500.

        z = NP.arange(zl,zr,zstep)
        ipoly_werte = [(x,0.,EzPoly(x, polyValues)) for x in z]

        ax  = plt.subplot(111)
        self.display(sf_data.Ez0_tab,'SFdata')
        self.display(ipoly_werte,'Polydata')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(input_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
        # TTF calculations
        v0  = V0(polyValues,zintval)
        t   = T(polyValues,k,zintval)
        s   = S(polyValues,k,zintval)
        tp  = Tp(polyValues,k,zintval)
        sp  = Sp(polyValues,k,zintval)
        DEBUG_OFF('V0 {}'.format(v0))
        DEBUG_OFF('T(k) {}'.format(t))
        DEBUG_OFF("T'(k) {}".format(tp))
        DEBUG_OFF('S(k) {}'.format(s))
        DEBUG_OFF("S'(k) {}".format(sp))                                                
    def test5(self):
        print("\b----------------------------------------test5")
        input_file='SF/PILL-2CM.TBL'
        input_file='SF/SF_WDK2g44.TBL'
        input_file='SF/CAV-FLAT-R135-L31.TBL'
        EzPeak  = 10
        IntgInterval = 0.
        sfdata = SFdata.field_data(input_file,EzPeak=EzPeak,IntgInterval=IntgInterval)
        print("peak:{} -- average:{} -- average/peak {}".format(sfdata.EzPeak,sfdata.EzAvg,sfdata.EzAvg/sfdata.EzPeak))
    def test6(self):
        print("\b----------------------------------------test6")
        input_file='SF/CAV-FLAT-R135-L31.TBL'
        IntgInterval = 0.
        EzPeak  = [0,1,5,10]
        wbool = True
        ix = 0
        while(wbool):
            print(f"\neffective gap for EzPeak {EzPeak[ix]}, IntgInterval CONST")
            sfdata     = SFdata.field_data(input_file,EzPeak=EzPeak[ix],IntgInterval=IntgInterval)
            EzAvg      = sfdata.EzAvg
            polyValues = sfdata.polyValues
            zEPoly(EzAvg,polyValues)
            ix += 1
            if ix > 3: wbool = False
        return
    def test7(self):
        print("\b----------------------------------------test7")
        input_file='SF/CAV-FLAT-R135-L31.TBL'
        IntgInterval = [0,3,9,18]
        EzPeak  = 0
        wbool = True
        ix = 0
        while(wbool):
            print(f"\neffective gap for EzPeak CONST, IntgInterval {IntgInterval[ix]}")
            sfdata     = SFdata.field_data(input_file,EzPeak=EzPeak,IntgInterval=IntgInterval[ix])
            EzAvg      = sfdata.EzAvg
            polyValues = sfdata.polyValues
            zEPoly(EzAvg,polyValues)
            ix += 1
            if ix > 3: wbool = False
        return

if __name__ == '__main__':
    # unittest.main()
    tests = TestEz0Methods()
    tests.test0()    
    tests.test1()
    tests.test2()
    tests.test3()
    tests.test4()
    tests.test5()
    tests.test6()
    tests.test7()
    