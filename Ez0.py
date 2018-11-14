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
import matplotlib.pyplot as plt
import numpy as NP
from math import sin,tan,pi,exp,fmod,cos
from collections import namedtuple

from setutil import PARAMS,DEBUG,Proton

# Polyval: 
# polynomial approximation for E(z,r=0), z in interval [zl,zr]: see (4.4.1) A.Shishlo/J.Holmes
Polyval = namedtuple('Polyval',['zl','z0','zr','dz','b','a','E0','coeff'])

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF
DEBUG_TEST2  = DEBUG_OFF
DEBUG_TEST3  = DEBUG_OFF

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

def EzAvg(poly):
    """ Average E-field in gap """
    sum = 0.
    N   = len(poly)
    for n in range(N):
        dz  = poly[n].dz      # [cm]
        v0  = V0n(poly,n)     # [MV]
        Eav = v0/(dz*1.e-2)   # [Mv/m]
        sum = sum + Eav
        pass
    Eav = sum/N               # EzAvg [MV/m]
    return Eav

class DynacSteps(object):
    """ Produce DYNAC's steps and parts values.
        DYNAC integrates the gap in serveral steps depending on the Superfish data.
        Each step is split into 4 parts. 
        Each part has a (z,time) tuple.
        The tuples are arranged in 2 matrices: one for z-values and one for time-values.
        Time values are calculated from velocity beta in units of lightspeed.
    """
    def __init__(self,gap,SFdata):
        self.zl = -gap*100./2.   # [cm]
        self.zr = -self.zl
        polyvals = []
        for poly in SFdata.Ez_poly:
            zil = poly.zl
            zir = poly.zr
            if zil < self.zl or zir > self.zr: 
                continue
            else:
                polyvals.append((poly.zl*1.e-2,poly.zr*1.e-2))    # all-in [m]
        polyvals = tuple(polyvals) 
        self.h = polyvals[0][1] - polyvals[0][0]
        z_parts = []
        for interval in polyvals:
            z0 = interval[0]
            z1 = z0 + self.h/4.
            z2 = z0 + self.h/2.
            z3 = z0 + (3*self.h)/4.
            z4 = z0 + self.h
            z_parts.append((z0,z1,z2,z3,z4))
        self.z_steps = tuple(z_parts)
        return
    
    def _z_step(self,num):
        return self.z_steps[num]
    
    def nsteps(self):
        return len(self.z_steps)

    def nparts(self):
        return 5
    
    def steplen(self):
        return self.h

    def _t_parts(self,beta,num):
        velocity = beta * 2.99792458e8
        t0 = self.z_steps[num][0]/velocity
        t1 = t0 + self.h/(4*velocity)
        t2 = t0 + self.h/(2*velocity)
        t3 = t0 +(3*self.h)/(4*velocity)
        t4 = t0 + self.h/velocity
        return (t0,t1,t2,t3,t4)
 
    def __call__(self,beta):
        nsteps = self.nsteps()
        nparts = self.nparts()
        zs = NP.zeros((nsteps,nparts))
        ts = NP.zeros((nsteps,nparts))
        for i in range(nsteps):
            z_step = self._z_step(i)
            t_step = self._t_parts(beta,i)
            for j in range(nparts):
                zs[i,j] = z_step[j]
                ts[i,j] = t_step[j]
        return zs,ts
            

class SFdata(object):
    ''' Cavity E(z,r=0) field profile: Superfish data  (normiert auf EzPeak)'''
    def __init__(self,input_file,EzPeak=1.):
        print('READING SF-DATA from "{}"'.format(input_file))
        self.input_file = input_file
        self.EzPeak     = EzPeak
        self.make_Ez_table()
        self.make_Ez_poly()
        self.EzAvg = EzAvg(self._poly)

    def make_Ez_table(self):
        """ read data and normalize """
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
        enorm = self.EzPeak/emax
        self._Ez0_tab = []
        for i in range(len(zp)):   # normalize and pack
            ep[i] = ep[i]*enorm
            dpoint = Dpoint(zp[i],rp[i],ep[i])
            self._Ez0_tab.append(dpoint)

    def make_Ez_poly(self):
        """ Calculate polynom coefficients from SF data """
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
                DEBUG_MODULE('make_Ez_poly:indexer(): (i,i+n,i+2*n)={},{},{}'.format(i,i+n,i+2*n))
                yield((i,i+n,i+2*n))
        
        self._poly = []
        anz = 10        # interpolate SF-data with 'anz' polynomials od 2nd order
        for (il,i0,ir) in indexer(anz,len(self.Ez_table)):
            DEBUG_MODULE('make_Ez_poly(): (il,i0,ir) ',((il,i0,ir)))
            zl = self.Ez_table[il].z
            z0 = self.Ez_table[i0].z
            zr = self.Ez_table[ir].z
            El = self.Ez_table[il].Ez
            E0 = self.Ez_table[i0].Ez
            Er = self.Ez_table[ir].Ez
            dz = z0-zl
            b  = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
            a  = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
            pval = Polyval(zl,z0,zr,dz,b,a,E0,0.)
            self._poly.append(pval)

    def Ez0t(self, z, t, omega, phi):
        """E(z,0,t): time dependent field value at location z"""
        res = Ipoly(z,self.Ez_poly) * cos(omega*t+phi)
        return res

    def dEz0tdt(self, z, t, omega, phi):
        """dE(z,0,t)/dt: time derivative of field value at location z"""
        res = - omega * Ipoly(z,self.Ez_poly) * sin(omega*t+phi)
        return res

    @property
    def Ez_table(self):
        """List(Dpoint) of SuperFish data points"""
        return self._Ez0_tab
    @property
    def Ez_poly(self):
        """List(Polyval) of polygon approximations"""
        return self._poly
        
def pre_plt(input_file):
    """ prepare plot """
    ax  = plt.subplot(111)
    ax.set_ylabel('Ez0 [MV/m]')
    ax.set_xlabel('z [cm]')
    ax.set_title(input_file)
    return ax

def post_plt(ax):
    """ finish plot """
    plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

def displayLR(table,legend):
    """ display left (L) and right (R) table data """
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
    """SuperFish (SF) data"""
    Ez0_tab = SFdata(input_file).Ez_table
    display(Ez0_tab,'SF-raw')
    
def test1():
    '''Gauss'che Normalverteilung (NG)'''
    print('----------------------------TEST1---')
    gap = 4.4
    z = NP.arange(0.,gap,gap/500.)
    sigma = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigma,0.)) for x in z]
    displayLR(Ez0_tab,'NG')

            # EXPERIMENTAL
    # EzZ_tab = [(x,0.,Z(x,sigma,0.01 )) for x in z]
    # display(EzZ_tab,'Z')
            # END-EXPERIMENTAL

            # Overlap?
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
            # END-Overlap?

def test2():
    '''Second order polynomial fit with 3 point formula to NG formula'''
    print('----------------------------TEST2---')
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

    z = NP.linspace(zl,zr,2*anz+1)
    Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
    # display(Ez0_tab,'slice')
    poly  = Kpoly(z,sigma,0.,E0*1.)

    zstep = (zr-zl)/500.
    z = NP.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez0_tab,'NG-poly')

            # Overlap?
    # poly1 = Kpoly(z,sigma,gap,E0*1.)  # next cavity
    # Ez1_tab = [(x,0.,Ipoly(x,poly1)) for x in z]
    # display(Ez1_tab,'poly1')
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # display(Ez2_tab,'poly2')
            # END-Overlap?

    # TTF calculations
    zintval = (zl,zr)
    v0  = V0(poly,zintval)
    t   = T(poly,k,zintval)
    s   = S(poly,k,zintval)
    tp  = Tp(poly,k,zintval)
    sp  = Sp(poly,k,zintval)
    DEBUG_TEST2('V0',v0)
    DEBUG_TEST2('T(k)',t)
    DEBUG_TEST2("T'(k)",tp)
    DEBUG_TEST2('S(k)',s)
    DEBUG_TEST2("S'(k)",sp)

def test3(input_file):
    '''Second order polynomial fit with 3 point formula to SF data'''
    print('----------------------------TEST3---')
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)*1.e-2    # [1/cm]
    gap = 4.8
    zl  = -gap/2.
    zr  = +gap/2.
    zintval = (zl,zr)

    gap_data = SFdata(input_file,EzPeak=1.)
    poly     = gap_data.Ez_poly
    display(gap_data.Ez_table,'SF-slice')

    zstep  = (zr-zl)/500.
    z      = NP.arange(zl,zr,zstep)
    Ez_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez_tab,'SF-poly')

    # TTF calculations
    v0  = V0(poly,zintval)
    t   = T(poly,k,zintval)
    s   = S(poly,k,zintval)
    tp  = Tp(poly,k,zintval)
    sp  = Sp(poly,k,zintval)
    DEBUG_TEST3('V0',v0)
    DEBUG_TEST3('T(k)',t)
    DEBUG_TEST3("T'(k)",tp)
    DEBUG_TEST3('S(k)',s)
    DEBUG_TEST3("S'(k)",sp)

def test4(input_file):
    print('----------------------------TEST4---')
    sf_data = SFdata(input_file)
    gap = 0.048
    beta = 0.1
    steps = DynacSteps(gap,sf_data)
    zs,ts = steps(beta)
    h = steps.steplen()
    # show
    print('gap [m]= ',gap,' beta= ',beta,'  h[m]= ',h)
    for ns in range(steps.nsteps()):
        for np in range(steps.nparts()):
            print('z [m]= ',zs[ns,np],'\tt[sec]= ',ts[ns,np])
        print()
    return
    
if __name__ == '__main__':
    input_file='SF_WDK2g44.TBL'
    ax = pre_plt(input_file)
    # test0(input_file)               # SF
    # test1()                         # NG
    # test2()                         # poly-fit to NG
    test3(input_file)               # poly-fit to SF
    post_plt(ax)
    test4(input_file)
