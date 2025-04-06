#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='11.0.2.4'
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
from math import exp,pi,tan,sin
from collections import namedtuple
from setutil import PARAMS,Proton,DEBUG_ON,DEBUG_OFF,wrapRED

# Polyval: polynomial approximation for E(z,r=0), z in interval [zl,zr]: see (4.4.1) A.Shishlo/J.Holmes
Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0'])
# Dpoint: Table data point - Ez0_tab is list(Dpoint)
Dpoint  = namedtuple('Dpoint',['z','R','Ez'])

def normGauss(x,sig,mu):
     # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res
def gaussPoly(z,sigma,mu,E):
    #  Polynom coefficients from Gauss Normalverteilung
    poly = []
    anz = len(z)
    for i in range(0,anz-2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  = z0-zl
        El  = E*normGauss(zl,sigma,mu)
        E0  = E*normGauss(z0,sigma,mu)
        Er  = E*normGauss(zr,sigma,mu)
        b = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation
        a = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
        pval = Polyval(zl,z0,zr,dz,b,a,E0)
        poly.append(pval)
    return poly
def EPoly(z,poly):
    """ Interpolation using the polynomial fit """
    ix = -1
    for i in range(len(poly)):
        zl = poly[i].zl
        zr = poly[i].zr
        if zl <= z and z <= zr:
            ix = i
            break
    if ix < 0:
        raise RuntimeError('EPoly(): arg out of range! {}'.format(z))

    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res
def zPoly(EzAvg,polyValues):
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
    e2 = EPoly(z2,polyValues)
    e1 = EPoly(z1,polyValues)
    e0 = EPoly(z0,polyValues)
    DEBUG_OFF(f"z0 {z0},z1 {z1},z2 {z2},e0 {e0},e1 {e1},e2 {e2},")
    #   use NP.polyfit(...,deg=2)
    pfit = NP.polyfit([e0,e1,e2],[z0,z1,z2],deg=2)
    zeff = pfit[0]*EzAvg**2 + pfit[1]*EzAvg + pfit[2]
    DEBUG_OFF(f"input EzAvg {EzAvg}, fitted EzAvg(zeff) {EPoly(zeff,polyValues)}, zeff {zeff}")
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
    (zl,zr) = zintval
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
    (zl,zr) = zintval
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
    (zl,zr) = zintval
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
    (zl,zr) = zintval
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
    (zl,zr) = zintval
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        if zil < zl or zir > zr: continue
        sp.append(Spn(poly,k,i))
    return sp

class SFdata(object):
    """ 
    Cavity E(z,r=0) a.k.a. Ez0 field profile from SF *.TBL file.
    Raw data can be scaled for max ordinate (EzPeak) and max abscissa (L)
    Data from SF are defined on positive interval z in [0,+L/2]
    Final data will be symmetric, i.e. E(-z)=E(+z) for z in [-L/2,+L/2]
    """
    def __init__(self,TBL_file):
        print('READING SF-DATA from "{}"'.format(TBL_file))
        # raw members start with underline (will never be modified)
        self._Ez0_tab  = None  # _Ez0_tab: unscaled table on [-L/2,+L/2]
        self._EzPeak   = None  # peak Ez
        self._L        = None  # definition interval [-L/2,+L/2] in [cm]
        self._EzAvg    = None  # average Ez
        # scaled data will go here
        self.Ez0_tab   = None
        self.EzPeak    = None
        self.L         = None  # definition interval [-L/2,+L/2] in [cm]
        self.EzAvg     = None
        self.polies    = None  # list(Polyvals)
        self.HE_Gap    = None  # equivalent hard edge gap
        self.HE_EzPeak = None  # equivalen hard edge Ez

        self.TBL_file = TBL_file
        self.readRawData()  

    def readRawData(self):
        """ read raw raw table, mirror data on ordinate and calc table's average (EzAvg_raw) """
        zp = []; rp = []; ep = []
        leading  =  41      # nbof leading lines to skip
        trailing  = 2       # nbof trailing lines to skip
        with open(self.TBL_file,'r') as f:
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

        # construct symmetric table for interval [-L/2,+L/2]
        # (E0,z0)=(ep[0],zp[0])
        # (ER,zL)=(ep[-1],zp[-1])
        EzPeak=El=ep[0]
        EzMin=Er=ep[-1]
        zl = zp[0]
        zr = zp[-1]
        L  =(zr-zl)*2
        dz = zp[1]-zp[0]

        zprev = [-x for x in reversed(zp[1:])]
        zp = zprev+zp
        rprev = [x for x in reversed(rp[1:])]
        rp = rprev+rp
        eprev = [x for x in reversed(ep[1:])]
        ep = eprev+ep
        N = len(zp)
        raw_tab = [Dpoint(zp[i],rp[i],ep[i]) for i in range(N)]
        DEBUG_OFF('raw_tab',raw_tab)
        _EzAvg = float((NP.trapezoid(ep,dx=dz)/L))
        DEBUG_OFF(f'[(Ez,z)]=[({EzPeak},{zl}),({EzMin},{zr})], L={L}, dz={dz}, EzAvg ={_EzAvg}')

        # raw data from SF will never be modified!
        self._Ez0_tab = raw_tab #NOTE: defined on  [-L/2,+L/2]
        self._L       = L       #NOTE: _L is length of [-L/2,+L/2]
        self._EzPeak  = EzPeak
        self._EzAvg   = _EzAvg
        # self.Ez0_tab  = raw_tab
        # self.L        = L
        # self.EzPeak   = E0
        return

    # class variables
    instances = {}    # dict of sigleton objects

    @classmethod
    def InstanciateAndScale(cls,TBL_file,EzPeak=0.,L=0.):
        """ 
        Return a singleton scaled SFdata object bound to TBL_file.
        IN: 
            EzPeak peak field [MV/m]
            L [cm]: interval z in [-L/2,+L/2] on which Ez0 is defined.
        OUT: 
            SFdata object 
        """
        # use a key made from <SFdata file name>%<Ezpeak value>%<LInterv value> to retrieve the corresponding field table instance
        instance_key = '{}%{}%{}'.format(TBL_file,EzPeak,L)
        # instance with instance_key exists? get it; else None
        instance = cls.instances.get(instance_key)
        if instance == None:
            instance = SFdata(TBL_file)   # new scaled instance
            if EzPeak != 0. and L != 0.:
                instance.scaleEzTable(EzPeak,L)
            elif EzPeak != 0. and L == 0.:
                instance.scaleEzTable(EzPeak,instance._L)
            elif EzPeak == 0. and L != 0.:
                instance.scaleEzTable(instance._EzPeak,L)
            elif EzPeak == 0. and L == 0.:
                instance.scaleEzTable(instance._EzPeak,instance._L)
            SFdata.instances[instance_key] = instance
            DEBUG_OFF(SFdata.instances)
            DEBUG_OFF('self.Ez0_tab',instance.Ez0_tab,'================= EOF Ez0_tab')
            instance.makeEPoly()
            DEBUG_OFF('makeEzPoly: {} poly intervals'.format(len(instance.polies)))
        return instance

    def Ez0t(self, z, t, omega, phis):
        """E(z,0,t): time dependent field value at location z"""
        res = EPoly(z,self.polyValues) * cos(omega*t+phis)
        return res
    def dEz0tdt(self, z, t, omega, phis):
        """dE(z,0,t)/dt: time derivative of field value at location z"""
        res = - omega * EPoly(z,self.polyValues) * sin(omega*t+phis)
        return res
    def scaleEzTable(self,EzPeak,L):   #NOTE: L is length [-L/2,+L/2] in [cm]
        self.EzPeak = EzPeak
        self.L = L
        Ez0_tab = []
        Ex = self._EzPeak
        Lx  = self._L
        for ix in range(len(self._Ez0_tab)):
            z = self._Ez0_tab[ix].z*(L/Lx)       # scale z-axis
            r = self._Ez0_tab[ix].R*(L/Lx)       # scale R same as z
            e = self._Ez0_tab[ix].Ez*(EzPeak/Ex) # sclae Ez-axis
            Ez0_tab.append(Dpoint(z,r,e))
        self.EzAvg = self._EzAvg*(EzPeak/Ex)
        self.Ez0_tab = Ez0_tab
        return
    def makeEPoly(self):
        # Polynomial fits to raw data according to Shislo&Holmes
        # In here M adjacent raw intervals are taken as one single interval and fitted with a polynomial
        sf_tab = self.Ez0_tab    # SF-table scales applied
        N      = len(sf_tab)     # stuetzpunkte in raw table
        M      = 8               # raw intervalls/poly interval (must be even number [2,4,6,8,....])
        polies = []              # polies: list(Polyvals)

        DEBUG_OFF('makeEzPoly: raw function values: {} in {}'.format(N,range(N-1)))
        DEBUG_OFF('makeEzPoly: first is sf_tab[{:3}]..{}'.format(0,sf_tab[0]))
        DEBUG_OFF('makeEzPoly: last is  sf_tab[{:3}]..{}'.format(N-1,sf_tab[N-1]))
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
            DEBUG_OFF('Ez0_poly::SFdata::makeEzPoly: (il,i0,ir) ({:3},{:3},{:3}),  (zl,z0,zr,E0) ({:6.3f},{:6.3f},{:6.3f},{:6.3f})'.format(il,i0,ir,zl,z0,zr,E0))
        DEBUG_OFF('makeEzPoly: {} poly intervals'.format(len(polies)))
        self.polies = polies
        return
    def hardEdge(self,gap):
        """ calculate equivalent hard edge parameters
        IN:
            gap [cm]: full hard edge gap size
        OUT:
            (HE_Gap,HE_EzPeak) [cm,MV/m]: (hard edge gap, hard edge peak field)
        """
        area = self.EzAvg*self.L
        self.HE_EzPeak = area/gap
        self.HE_Gap = gap
        return (self.HE_Gap,self.HE_EzPeak)

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
    def display(cls,table,legend,linestyle="solid"):
        z   = [+float(x[0]) for x in table]
        Ez  = [+float(x[2]) for x in table]
        plt.plot(z,Ez,label=legend,linestyle=linestyle)
    def test0(self):
        print("\b---------------------------------test0 (singleton)")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        (EzPeak,L) = (0.0,0.0)
        sfdata1 = SFdata.InstanciateAndScale(TBL_file,*(EzPeak,L))
        print('   ',(EzPeak,L),sfdata1.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,L) = (4.0,0.0)
        sfdata3 = SFdata.InstanciateAndScale(TBL_file,*(EzPeak,L))
        print('   ',(EzPeak,L),sfdata3.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,L) = (0.,4.)
        sfdata4 = SFdata.InstanciateAndScale(TBL_file,*(EzPeak,L))
        print('   ',(EzPeak,L),sfdata4.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        # TBL_file='SF/SF_WDK2g22.TBL'
        TBL_file='SF/ALCELI-750.0-2.02.TBL'
        sfdata5 = SFdata.InstanciateAndScale(TBL_file)
        print('   ',(EzPeak,L),sfdata5.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        (EzPeak,L) = (4.0,4.0)
        sfdata6 = SFdata.InstanciateAndScale(TBL_file,*(EzPeak,L))
        print('   ',(EzPeak,L),sfdata6.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])

        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/ALCELI-750.0-2.02.TBL'
        (EzPeak,L) = (0.,0.)
        sfdata7 = SFdata.InstanciateAndScale(TBL_file)
        print('   ',(EzPeak,L),sfdata7.TBL_file,end="")
        print(' => SFdata.instances ',[id(x) for x in SFdata.instances])
    def test1(self):
        print("\b----------------------------------------test1")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        sfdata  = SFdata.InstanciateAndScale(TBL_file,EzPeak=1.5,L=5.)
        Ez0_tab        = sfdata._Ez0_tab    # raw
        Ez1_tab        = sfdata.Ez0_tab     # scaled
        sfdata.scaleEzTable(EzPeak=2.5,L=2.)
        Ez2_tab        = sfdata.Ez0_tab     # scaled

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'SF-raw')
        self.display(Ez1_tab,'SF-scaled.1')
        self.display(Ez2_tab,'SF-scaled.2')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(TBL_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test2(self):
        print("\b----------------------------------------test2")
        L = 4.4
        z = NP.arange(0.,L,L/500.)
        sigma = 1.14
        Ez0_tab = [(x,0.,normGauss(x,sigma,0.)) for x in z]

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
        L     = 4.8          # [cm] full interval length
        zl    = -L/2.        #left  interval boundary
        zr    = L/2.         #right interval boundary
        sigma = L/2./1.89    # sigma of normGauss (best fit with SF)
        # sigma = L/2./2.2   # sigma of normGauss (best fit with SF)
        E0    = 1.           # top of normGauss   (best fit with SF)

        z = NP.linspace(zl,zr,2*anz+1)
        Ez0_tab = [(x,0.,E0*normGauss(x,sigma,0.)) for x in z]
        # display(Ez0_tab,'slice')
        poly = gaussPoly(z,sigma,0.,E0*1.)

        zstep = (zr-zl)/500.
        z = NP.arange(zl,zr,zstep)
        Ez0_tab = [(x,0.,EPoly(x, poly)) for x in z]

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'NG-poly')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
    def test4(self):
        print("\b----------------------------------------test4")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        EzPeak    = 1.2
        sf_data   = SFdata.InstanciateAndScale(TBL_file,EzPeak=EzPeak,L=0.)
        polyValues= sf_data.polies
        particle  = Proton(tkin=100.)
        beta      = particle.beta
        c         = PARAMS['clight']
        freq      = 800.e6
        k         = 2*pi*freq/(c*beta)*1.e-2    # [1/cm]
        zr        = sf_data.L/2
        zl        = -zr
        zintval   = (zl,zr)
        zstep     = (zr-zl)/500.

        z = NP.arange(zl,zr,zstep)
        ipoly_werte = [(x,0.,EPoly(x, polyValues)) for x in z]

        ax  = plt.subplot(111)
        self.display(ipoly_werte,    legend='Polydata')
        self.display(sf_data.Ez0_tab,legend='SFdata',linestyle='dotted')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(TBL_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()
        # TTF calculations
        v0  = V0(polyValues,    zintval)
        t   = T( polyValues, k, zintval)
        s   = S( polyValues, k, zintval)
        tp  = Tp( polyValues,k, zintval)
        sp  = Sp( polyValues,k, zintval)
        DEBUG_OFF('V0',['{:.4g}'.format(x) for x in v0])
        DEBUG_OFF('T(k)',['{:.4g}'.format(x) for x in t])
        DEBUG_OFF('S(k)',['{:.4g}'.format(x) for x in s])
        DEBUG_OFF("T'(k)",tp)
        DEBUG_OFF("S'(k)",sp)                                         
    def test5(self):
        print("\b----------------------------------------test5")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        EzPeak  = 10
        L = 0.
        sfdata = SFdata.InstanciateAndScale(TBL_file,EzPeak=EzPeak,L=L)
        av2peak1 = sfdata.EzAvg/sfdata.EzPeak
        DEBUG_OFF("peak:{:.3f} -- average:{:.3f} -- average/peak {:.3f}".format(sfdata.EzPeak,sfdata.EzAvg,av2peak1))
        EzPeak = 4.5
        sfdata = SFdata.InstanciateAndScale(TBL_file,EzPeak=EzPeak,L=L)
        av2peak2 = sfdata.EzAvg/sfdata.EzPeak
        DEBUG_OFF("peak:{:.3f} -- average:{:.3f} -- average/peak {:.3f}".format(sfdata.EzPeak,sfdata.EzAvg,av2peak2))
        self.assertAlmostEqual(av2peak1,av2peak2,msg='average/peak',places=3)
    def test6(self):
        print("\b----------------------------------------test6")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        L = 0.
        EzPeak  = [0,1,5,15]
        zeff = []
        wbool = True
        ix = 0
        while(wbool):
            sfdata = SFdata.InstanciateAndScale(TBL_file,EzPeak=EzPeak[ix],L=L)
            print(f"effective gap for EzPeak {EzPeak[ix]}, L CONST")
            EzAvg      = sfdata.EzAvg
            polyValues = sfdata.polies
            zeff.append(zPoly(EzAvg,polyValues))
            if ix > 0:
                self.assertAlmostEqual(zeff[0],zeff[ix],msg='zeff',places=3)
            ix += 1
            if ix > 3: wbool = False
        
        return
    def test7(self):
        print("\b----------------------------------------test7")
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        L = [0,3,9,18]
        EzAvg = []
        zeff = []
        wbool = True
        ix = 0
        while(wbool):
            sfdata = SFdata.InstanciateAndScale(TBL_file,EzPeak=0,L=L[ix])
            print(f"effective gap for EzPeak CONST, L {L[ix]}")
            EzAvg.append(sfdata.EzAvg)
            polyValues = sfdata.polies
            zeff.append(zPoly(EzAvg[ix],polyValues))
            DEBUG_OFF(f"input EzAvg {EzAvg[ix]}, fitted EzAvg(zeff) {EPoly(zeff[ix],polyValues)}, zeff {zeff[ix]}")
            if ix > 0:
                self.assertAlmostEqual(EzAvg[ix-1],EzAvg[ix],msg='EzAvg',places=3)
                self.assertNotEqual(zeff[ix-1],zeff[ix],msg='zeff')
            ix += 1
            if ix > 3: wbool = False
        return
    def test8(self):
        print("\b----------------------------------------test8")
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        TBL_file='SF/SF_WDK2g22.TBL'
        TBL_file='SF/ALCELI-750.0-2.02.TBL'
        sfdata = SFdata.InstanciateAndScale(TBL_file)    # raw, no scaling
        reduc = 0.57   # reduce L by 57%
        gap = sfdata.L * reduc
        (HE_Gap, HE_EzPeak) = sfdata.hardEdge(gap)
        DEBUG_OFF(f'raw: (EzPeak,L)= ({sfdata.EzPeak:.3f},{sfdata.L:.3f}), (EzAvG,L)= ({sfdata.EzAvg:.3f},{sfdata.L:.3f}); hard edge: (He_EzPeak,HE_Gap)= ({HE_EzPeak:.3f},{HE_Gap:.3f})')
        self.assertAlmostEqual(HE_EzPeak*reduc,sfdata.EzAvg,msg='hard edge field',places=3)

        sfdata.scaleEzTable(1.,4.4)   # apply scaling
        gap = sfdata.L * reduc 
        (HE_Gap, HE_EzPeak) = sfdata.hardEdge(gap)
        DEBUG_OFF(f'raw: (EzPeak,L)= ({sfdata.EzPeak:.3f},{sfdata.L:.3f}), (EzAvG,L)= ({sfdata.EzAvg:.3f},{sfdata.L:.3f}); hard edge: (He_EzPeak,HE_Gap)= ({HE_EzPeak:.3f},{HE_Gap:.3f})')
        self.assertAlmostEqual(HE_EzPeak*reduc,sfdata.EzAvg,msg='hard edge field',places=3)

        sfdata.scaleEzTable(1.,4.4)  # same scaling, same instance
        reduc = 1. # reduce L by 0%
        gap = sfdata.L * reduc
        (HE_Gap, HE_EzPeak) = sfdata.hardEdge(gap)
        DEBUG_OFF(f'raw: (EzPeak,L)= ({sfdata.EzPeak:.3f},{sfdata.L:.3f}), (EzAvG,L)= ({sfdata.EzAvg:.3f},{sfdata.L:.3f}); hard edge: (He_EzPeak,HE_Gap)= ({HE_EzPeak:.3f},{HE_Gap:.3f})')
        self.assertAlmostEqual(HE_EzPeak*reduc,sfdata.EzAvg,msg='hard edge field',places=3)
    def test9(self):
        print("\b----------------------------------------test9")
        # print('test TTF calculation according to T.Wangler pp. 39 (2.14)')
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/SF_WDK2g22.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        sfdata = SFdata.InstanciateAndScale(TBL_file)    # raw, no scaling
        Ez0_tab = sfdata.Ez0_tab
        DEBUG_ON('Ez0 table',Ez0_tab)
        pass
    def test9(self):
        print("\b----------------------------------------test9")
        print(wrapRED('Gap-Scaling as used in TT28_ttf.yml'))
        # TBL_file='SF/PILL-2CM.TBL'
        # TBL_file='SF/SF_WDK2g44.TBL'
        TBL_file='SF/CAV-FLAT-R135-L32.TBL'
        # TBL_file='SF/ALCELI-750.0-2.02.TBL'
        sfdata  = SFdata.InstanciateAndScale(TBL_file,EzPeak=0,L=0.0942)
        Ez0_tab        = sfdata._Ez0_tab    # raw
        Ez1_tab        = sfdata.Ez0_tab     # scaled

        ax  = plt.subplot(111)
        self.display(Ez0_tab,'SF-raw')
        self.display(Ez1_tab,'SF-scaled.1')
        ax.set_ylabel('Ez0 [MV/m]')
        ax.set_xlabel('z [cm]')
        ax.set_title(TBL_file)
        plt.legend(loc='upper right',fontsize='x-small')
        plt.show()



if __name__ == '__main__': 
    # unittest.main()
    tests = TestEz0Methods()
    # tests.test0()    
    # tests.test1()
    # tests.test2()
    # tests.test3()
    # tests.test4()
    # tests.test5()
    # tests.test6()
    # tests.test7()
    # tests.test8()
    tests.test9()
    