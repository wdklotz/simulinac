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
from math import sin,cos,tan,radians,degrees,sqrt
from math import pi as PI
from copy import copy
import numpy as NP

from setutil import PARAMS,DEBUG,I0,I1,tblprnt,arrprnt
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from Ez0 import SFdata

# DEBUG__*
def DEBUG_ON(string,arg='',end='\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg='',end='\n'):
    pass
DEBUG_TEST0    = DEBUG_ON
DEBUG_TEST1    = DEBUG_ON
DEBUG_SLICE    = DEBUG_OFF
DEBUG_TTF_G    = DEBUG_OFF

twopi          = 2*PI

class _TTF_G(object):
    """Transition Time Factors RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247)"""
    def __init__(self, parent):
        def make_slices(parent, gap, SFdata, particle):
            """Slice the RF gap"""
            slices = []
            zl = -gap/2.*100.   # [m] --> [cm]
            zr = -zl
            E0z = 0.
            z = 0.
            for poly in SFdata.Ez_poly:
                zil = poly.zl
                zir = poly.zr
                if zil < zl or zir > zr: continue
                # instanciate _TTF_Gslices
                slice = _TTF_Gslice(parent, poly, particle)
                slices.append(slice)
            return slices

        def configure_slices(slices, phis, tkin):
            """adjust energy of slices"""
            next_phase = phis
            next_tkin  = tkin
            # Tklist = []         # helper to keep values to calculate min and max
            for slice in slices:
                # setting phase @ slice entrance
                slice.phis = next_phase
                # setting energy @ slice entrance
                slice.adjust_slice_parameters(next_tkin)
                # Tklist.append(slice.Tk)
                next_phase = slice.PHOUT       # slice OUT as next slice IN
                next_tkin  = slice.WOUT        # slice OUT as next slice IN
                DEBUG_TTF_G('_TTF_G: {}\n'.format(self),self.__dict__)
            deltaW  = next_tkin-tkin                        # total energy kick as sum over slices
            return deltaW

        # _TTF_G
        self.EzAvg    = parent.EzAvg
        self.gap      = parent.gap
        self.E0L      = self.EzAvg*self.gap
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.dWf      = parent.dWf
        self.lamb     = parent.lamb
        self.SFdata   = parent.SFdata
        self.matrix   = parent.matrix
        self.particle = parent.particle
        self.position = parent.position
        self.tkin     = self.particle.tkin
        if parent.SFdata == None:
            raise RuntimeError('_TTF_G: missing E(z) table - STOP')
            sys.exit(1)
        else:
             # slice the gap
            self.slices = make_slices(self, self.gap, self.SFdata, self.particle)
            # slice energy dependence
            self._deltaW = configure_slices(self.slices, self.phis, self.tkin)
            self._ttf = self._deltaW/(self.E0L*cos(self.phis)) if self.dWf == 1 else 1.
            # UPDATE linear NODE matrix with deltaW
            self.matrix[EKOO,DEKOO] = self._deltaW
            self._particlef = copy(self.particle)(self.particle.tkin + self._deltaW)
            # for test0()
            if DEBUG_TEST0 == DEBUG_ON:  parent['slices'] = self.slices

    # delegated parent properties
    @property
    def ttf(self):
        return self._ttf
    @property
    def deltaW(self):
        return self._deltaW
    @property
    def particlef(self):
        return self._particlef

    def map(self, i_track):
        """ Mapping from position (i) to (f )"""
        f_track = copy(i_track)
        # full map through sliced TTF-gap
        f_track = self._full_gap_map(self.slices, f_track)
        f_track[EKOO] += self._deltaW
        DEBUG_OFF('ttf-map ',f_track)
        return f_track

    def soll_map(self, i_track):
        si,sm,sf = self.position
        f_track = copy(i_track)
        f_track[EKOO] += self._deltaW
        f_track[SKOO]  = sm
        DEBUG_OFF('ttf-soll ',f_track)
        return f_track
        
    def _full_gap_map(self, slices, i_track):
        """ The wrapper to slice mappings """
        f_track = copy(i_track)
        for slice in slices:
            # map each slice with TTF 3-point gap-model
            f_track = slice.slice_map(f_track)
            # relativistic scaling. Is it needed?
            # z = f_track[ZKOO]
            # betai = self.particle.beta
            # tkin  = self.particle.tkin
            # betaf = Proton(tkin=tkin+self.deltaW).beta
            # z = betaf/betai*z
            # f_track[ZKOO] = z
        return f_track
        
class _TTF_Gslice(object):
    """ PyOrbit's Transit Time Factor RF-Gap Model """
    def __init__(self, parent, poly, particle):
        self.parent     = parent # the gap this slice belongs to
        self.freq       = parent.freq
        self.lamb       = parent.lamb
        self.particle   = copy(particle) # incoming SOLL particle
        self.poly       = poly # polynom interval: ACHTUNG: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.V0         = self._V0(self.poly)
        self.beta       = self.particle.beta
        self.gamma      = self.particle.gamma
        self.gb         = self.particle.gamma_beta
        self.k          = twopi/(self.lamb*self.beta)
        self.Tk         = self._T (self.poly, self.k)
        self.Tkp        = self._Tp(self.poly, self.k)
        self.phis      = None  # initialized in configure_slices
        self.WIN       = None  # initialized in adjust_slice_parameters
        self.WOUT      = None  # initialized in adjust_slice_parameters
        self.deltaW    = None  # initialized in adjust_slice_parameters
        self.PHIN      = None  # initialized in adjust_slice_parameters
        self.PHOUT     = None  # initialized in adjust_slice_parameters
        self.deltaPHI  = None  # initialized in adjust_slice_parameters

    def _T(self, poly, k):    # A.Shishlo/J.Holmes (4.4.6)
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        # DEBUG_SLICE('_TTF_Gslice:_T: (T,k)',(t,k))
        return t

    def _Tp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp

    def _V0(self, poly):    # A.Shishlo/J.Holmes (4.4.3)
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)*1.e-2       # [cm] --> [m]
        v0 = v0*E0*self.parent.dWf
        return v0

    def adjust_slice_parameters(self, tkin):
        """ Adjust energy-dpendent parameters for this slice """
        self.particle(tkin)
        self.beta    = self.particle.beta
        self.gamma   = self.particle.gamma
        self.gb      = self.particle.gamma_beta
        self.k       = twopi/(self.lamb*self.beta)
        self.Tk      = self._T (self.poly,self.k)
        self.Tkp     = self._Tp(self.poly,self.k)

        m0c2  = self.particle.e0
        m0c3  = self.particle.m0c3
        omeg  = twopi*self.freq
        i0    = 1.
        i1    = 0.
        WIN   = self.particle.tkin
        PHIN  = self.phis
        DW    = self.wout_minus_win(self.V0,i0,self.Tk,0.,PHIN)
        WOUT  = WIN+DW
        DPHI  = self.phiout_minus_phiin(self.V0*omeg/(m0c3*self.gb**3), self.gamma,0.,i0,i1,self.Tk,0.,self.Tkp,0.,PHIN)
        PHOUT = PHIN+DPHI

        self.WIN      = WIN
        self.WOUT     = WOUT
        self.deltaW   = DW
        self.PHIN     = PHIN
        self.PHOUT    = PHOUT
        self.deltaPHI = DPHI
        DEBUG_SLICE('_TTF_Gslice: {}\n'.format(self),self.__dict__)
        return 

    def wout_minus_win(self, fac, i0, tk, sk, phi):
        """Formel 4.3.1 A.Shishlo/J.Holmes"""
        return fac*i0*(tk*cos(phi)-sk*sin(phi))

    def phiout_minus_phiin(self, fac, gamma, r, i0, i1, tk, sk, tkp, skp, phi):
        """Formel 4.3.2 A.Shishlo/J.Holmes"""
        return  fac*i0*(tkp*sin(phi)+skp*cos(phi)+gamma*r*i1*(tk*sin(phi)+sk*cos(phi)))

    def slice_map(self, i_track):
        """Map through this slice from position (i) to (f)"""
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z~(phi-phis)
        zp       = i_track[ZPKOO]      # [5] dp/p~dT
        T        = i_track[EKOO]       # [6] kinetic energy soll
        S        = i_track[SKOO]       # [8] position soll

        c          = PARAMS['lichtgeschwindigkeit']
        m0c2       = self.particle.e0
        m0c3       = self.particle.m0c3
        omeg       = twopi*self.freq

        # energy parameters from SOLL
        betai      = self.particle.beta
        gammai     = self.particle.gamma
        gbi        = self.particle.gamma_beta

        pin        = -z*omeg/(c*betai) + self.PHIN            # phase  (i)  ~ (-z)
        win        = (zp*(gammai+1.)/gammai+1.)*self.WIN      # energy (i)  ~ (z')

        # energy parameters from particle ??
        # particle   = copy(self.particle)(tkin=win)  
        # betai      = particle.beta
        # gammai     = particle.gamma
        # gbi        = particle.gamma_beta

        k          = omeg/(c*betai)
        Tk         = self._T(self.poly,k)
        Tkp        = self._Tp(self.poly,k)
        r          = sqrt(x**2+y**2)            # radial coordinate
        K          = omeg/(c*gbi) * r
        i0         = I0(K)                      # bessel function I0
        i1         = I1(K)                      # bessel function I1
        fact       = self.V0*omeg/(m0c3*gbi**3)

        womwi = self.wout_minus_win(self.V0,i0,Tk,0.,pin) 
        wout  = win + womwi                              # energy (f)

        pompi = self.phiout_minus_phiin(fact,gammai,r,i0,i1,Tk,0.,Tkp,0.,pin)
        pout  = pin + pompi                              # phase (f)
        
        dp = +(pout-self.PHOUT)    # delta phase  (f)
        dw = +(wout-self.WOUT)     # delta energy (f)
        
        zf  = -dp*(c*betai)/omeg
        zpf = gammai/(gammai+1)*dw/self.WOUT

        f_particle = copy(self.particle)(tkin=self.WOUT)  # energy parameters from SOLL
        # f_particle = copy(self.particle)(tkin=wout)     # energy parameters from PARTICLE ??
        gbf        = f_particle.gamma_beta

        fact = self.V0/(m0c2*gbi*gbf)*i1
        if r > 0.:
            xp = gbi/gbf*xp-x/r*fact*Tk*sin(pin)
            yp = gbi/gbf*yp-y/r*fact*Tk*sin(pin)
        elif r == 0.:
            xp = gbi/gbf*xp
            yp = gbi/gbf*yp

        f_track = NP.array([x,xp,y,yp,zf,zpf,T,1.,S,1.])
        return f_track

def test0():
    import elements as ELM
    from bunch import Tpoint, Track
    
    print('-----------------------------------TEST 0----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']*1.8055 # [Mv/m] EzPeak/EzAvg fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,EzPeak)
    
    ttfg = ELM.RFG(gap=0.048,SFdata=SF_tab,mapping='ttf')
    tkin = 50.
    ttfg.adjust_energy(tkin=tkin)
    if DEBUG_TEST0 == DEBUG_ON:
        print('TTFG: ttfg.__dict__',ttfg.__dict__)      # for DEBUGGING
        slices = ttfg['slices']
        for slice in slices:
            print('_TTF_Gslice: slice\n',slice.__dict__)      # for DEBUGGING
            pass
    else:
        pass

    z = 1.e-3
    x=y=1.e-2
    T = tkin
    # track-point fields:              x   x'  y  y'  z   z'  T  1   S   1
    tpoint = Tpoint(point = NP.array([ x,  0., y, 0., z,  0., T, 1., 0., 1.]))
    track = Track()
    track.addpoint(tpoint)
    ti = track.getpoints()[-1]
    for i in range(1):
        DEBUG_TEST0('MAP:\n',track.getpoints()[-1].as_str())
        tf = ttfg.map(ti())
        tpf = Tpoint(tf)
        track.addpoint(tpf)
        DEBUG_TEST0('MAP:\n',track.getpoints()[-1].as_str())
        ttfg.adjust_energy(tf[EKOO])    #enery adaptation
        ti = tpf

def test1():
    import elements as ELM
    from bunch import Tpoint, Track
    
    print('-----------------------------------TEST 1----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']
    SF_tab = SFdata(input_file,EzPeak)
    
    ttfg = ELM.RFG(gap=0.048,SFdata=SF_tab,mapping='ttf')
    tkin = 150.
    ttfg.adjust_energy(tkin=tkin)
    DEBUG_TEST1('TTFG: ttfg.__dict__',ttfg.__dict__)      # for DEBUGGING
    
    von = 0.
    bis = 1.
    anz = 20
    delta = (bis-von)/anz
    x=xp=y=yp=z=zp=0.0
    x = 1.e-1
    z = von
    T = 0.
    # start                           x   x'  y  y'  z   z'  T  1   S   1
    start = Tpoint(point = NP.array([ x,  xp, y, yp, z,  zp, T, 1., 0., 1.]))
    track = Track()
    track.addpoint(start)
    ti = track.getpoints()[-1]
    for i in range(anz+1):
        tf = ttfg.map(ti())
        tpf = Tpoint(tf)
        track.addpoint(tpf)
        z += delta
        tf[4] = z
        ti = Tpoint(tf)
    DEBUG_TEST1('TRACK-POINTS:\n',track.as_table())

if __name__ == '__main__':
    test0()
    test1()