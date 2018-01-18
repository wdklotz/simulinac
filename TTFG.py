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
import sys
from math import sin,cos,tan,radians,degrees,sqrt
from math import pi as PI
from copy import copy
import numpy as NP

from setutil import PARAMS,DEBUG,I0,I1,tblprnt,arrprnt
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from Ez0 import SFdata

def DEBUG_ON(string,arg='',end='\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg='',end='\n'):
    pass
## DEBUG__*
DEBUG_TEST0    = DEBUG_ON
DEBUG_TEST1    = DEBUG_ON
DEBUG_MAP      = DEBUG_OFF
DEBUG_SOLL_MAP = DEBUG_OFF
DEBUG_SLICE    = DEBUG_OFF

twopi          = 2*PI

class _TTF_G(object):
    """Transition Time Factors RF Gap-Model (A.Shishlo ORNL/TM-2015/247)"""
    def __init__(self, parent):
        def _make_slices():
            """Slice the RF gap"""
            slices = []
            zl = -self.gap/2.*100.   # [m] --> [cm]
            zr = -zl
            E0z = 0.
            z = 0.
            for poly in self.SFdata.Ez_poly:
                zil = poly.zl
                zir = poly.zr
                if zil < zl or zir > zr: continue
                slice = _TTF_Gslice(self,poly,self.particle)  # instanciate _TTF_Gslices
                slices.append(slice)
            return slices

        def adjust_slice_energy(phis, tkin):
            """adjust energy of slices"""
            next_phase = phis
            next_tkin  = tkin
            Tklist = []         # helper to keep values to calculate min and max
            for slice in self.slices:
                slice.setSollPhase(next_phase) # NOTE!: phase  @ slice entrance set here
                slice.adjust_slice_parameters(next_tkin) # NOTE!: energy @ slice entrance set here
                Tklist.append(slice.Tk)
                next_phase = slice.PHOUT       # NOTE!: slice OUT becomes next slice IN
                next_tkin  = slice.WOUT        # NOTE!: slice OUT becomes next slice IN
                DEBUG_SLICE('SLICE: {}\n'.format(self),self.__dict__)
            deltaW   = next_tkin-tkin          # total energy kick of this gap
#todo:  @@@@ TODO should calc the average @@@@
            self.tr  = min(Tklist)
            return deltaW

        self.u0       = parent.u0
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.label    = parent.label
        self.gap      = parent.gap
        self.mapping  = parent.mapping
        self.dWf      = parent.dWf
        self.lamb     = parent.lamb
        self.SFdata   = parent.SFdata
        self.matrix   = parent.matrix
        self.particle = parent.particle
        self.position = parent.position
        if self.SFdata == None:
            raise RuntimeError('_TTF_G: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.slices = _make_slices()             # slice the gap
            self.deltaW = adjust_slice_energy(self.phis, self.particle.tkin)
            self.matrix[EKOO,DEKOO] = self.deltaW    # nicht vwergessen! set my deltaW in linear map
            if DEBUG_TEST0 == DEBUG_ON:  parent['slices'] = self.slices   # for test0()

    def _full_gap_map(self, i_track):  # the wrapper to slice mappings
        self.dbTab1Rows  = []          # for DEBUGGING
        self.dbTab1Headr = []          # for DEBUGGING
        self.dbTab2Rows  = []          # for DEBUGGING
        self.dbTab2Headr = []          # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:      # for DEBUGGING
            self.dbTab1Headr = ['pout','pin','pout-pin','dp=pout-POUT','wout','win','wout-win','WOUT','dw=wout-WOUT','qV0*10^3']
            self.dbTab2Headr = ['x*10^3','xp*10^3','y*10^3','yp*10^3','z*10^3','zp*10^3','r*10^3','Tk',"Tkp",'i0-1','i1']

        for cnt,slice in enumerate(self.slices):
            # DEBUG_MAP('MAP: ttfg-map: {} tkin {} '.format(cnt,self.particle.tkin),slice)
            f_track = slice._slice_map(i_track)    # map slice with TTF 3-point gap-model
            i_track = f_track

                # relativistic scaling. Is it needed?
                # z = f_track[ZKOO]
                # betai = self.particle.beta
                # tkin  = self.particle.tkin
                # betaf = Proton(tkin=tkin+self.deltaW).beta
                # z = betaf/betai*z
                # f_track[ZKOO] = z

        DEBUG_MAP('MAP: ttfg-map: track through slices:\n',tblprnt(self.dbTab1Headr,self.dbTab1Rows))
        DEBUG_MAP('MAP: ttfg-map: track through slices:\n',tblprnt(self.dbTab2Headr,self.dbTab2Rows))

        # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:
            f = f_track.copy()
            for i in range(len(f_track)-4):
                f[i] =f[i]*1.e3
            arrprnt(f,fmt='{:6.3g},',txt='ttf_map: ')

        return f_track

    def map(self, i_track):
        """Mapping from position (i) to (f)"""
        f_track = self._full_gap_map(i_track)   # NOTE: use local map with sliced TTFGap
        self.particlef = copy(self.particle)(f_track[EKOO])

        # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:
            f = f_track.copy()
            for i in range(len(f_track)-4):
                f[i] =f[i]*1.e3
            arrprnt(f,fmt='{:6.3g},',txt='ttf_map: ')
        return f_track

        # Transition Time Factors RF Gap Model

class _TTF_Gslice(object):
    def __init__(self, parent, polyval, particle):
        self.parent     = parent           # the element this slice is part off
        self.freq       = parent.freq
        self.particle   = copy(particle)   # incoming soll particle
        self.polyval    = polyval # polynom interval: ACHTUNG: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.V0         = self._V0(self.polyval)
        self.beta       = self.particle.beta
        self.gamma      = self.particle.gamma
        self.gb         = self.particle.gamma_beta
        self.k          = twopi/(self.parent.lamb*self.beta)
        self.Tk         = self._T (self.polyval,self.k)
        self.Tkp        = self._Tp(self.polyval,self.k)
        self.phis       = None  # initialize before using!
        self.WIN        = None  # initialize before using!
        self.WOUT       = None  # initialize before using!
        self.deltaW     = None  # initialize before using!
        self.PHIN       = None  # initialize before using!
        self.PHOUT      = None  # initialize before using!
        self.deltaPHI   = None  # initialize before using!

    def setSollPhase(self, value): # must be called after instantiation and before usage
        self.phis = value

    def _T(self, poly, k):    # A.Shishlo (4.4.6)
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        # DEBUG_MAP('MAP: (T,k)',(t,k))
        return t

    def _Tp(self, poly, k):   # A.Shishlo (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp

    def _V0(self, poly):    # A.Shishlo (4.4.3)
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)*1.e-2       # [cm] --> [m]
        v0 = v0*E0*self.parent.dWf
        return v0

    def mapSoll(self, i_track):
        """Mapping of soll track through a slice from position (i) to (f)"""
        DEBUG_SOLL_MAP('mapSoll',NP.array2string(i_track))
        f_track = copy(i_track)   # make a copy to prevent i_track to be overwritten
        f_track[EKOO] += self.deltaW
        return f_track

    def adjust_slice_parameters(self, tkin):
        """Adjust energy and -dpendent parameters for this slice"""
        self.particle(tkin)
        self.beta    = self.particle.beta
        self.gamma   = self.particle.gamma
        self.gb      = self.particle.gamma_beta
        self.k       = twopi/(self.parent.lamb*self.beta)
        self.Tk      = self._T (self.polyval,self.k)
        self.Tkp     = self._Tp(self.polyval,self.k)

        c     = PARAMS['lichtgeschwindigkeit']
        m0c3  = self.particle.e0*c
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
        DEBUG_SLICE('SLICE {}\n'.format(self),self.__dict__)
        return 

    def wout_minus_win(self, conv, i0, tk, sk, phi):
        """Formel 4.3.1 A.Shishlo"""
        return conv*i0*(tk*cos(phi)-sk*sin(phi))

    def phiout_minus_phiin(self, fac, gamma, r, i0, i1, tk, sk, tkp, skp, phi):
        """Formel 4.3.2 A.Shishlo"""
        return  fac*i0*(tkp*sin(phi)+skp*cos(phi)+gamma*r*i1*(tk*sin(phi)+sk*cos(phi)))

    def _slice_map(self, i_track):
        """Map through this slice from position (i) to (f)"""
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z~(phi-phis)
        zp       = i_track[ZPKOO]      # [5] dp/p~dT
        T        = i_track[EKOO]       # [6] summe aller dT
        s        = i_track[SKOO]       # [8] summe aller laengen

        c          = PARAMS['lichtgeschwindigkeit']
        m0c2       = self.particle.e0
        m0c3       = m0c2*c
        omeg       = twopi*self.freq

        be         = self.particle.beta
        ga         = self.particle.gamma
        pin        = -z*omeg/(c*be) + self.PHIN       # phase @ (i)  ~ (-z)
        win        = (zp*(ga+1.)/ga+1.)*self.WIN      # energy @ (i) ~ (zp)

        # particle   = copy(self.particle)(tkin=win)  # energy parameters from particle
        particle   = self.particle                    # energy parameters from soll
        be         = particle.beta
        ga         = particle.gamma
        gb         = particle.gamma_beta
        k          = omeg/(c*be)
        Tk         = self._T (self.polyval,k)
        Tkp        = self._Tp(self.polyval,k)
        r          = sqrt(x**2+y**2)                # radial coordinate
        Kbess      = omeg/(c*gb) * r
        i0         = I0(Kbess)                      # bessel function
        i1         = I1(Kbess)                      # bessel function
        fact       = self.V0*omeg/(m0c3*gb**3)

        wowi = self.wout_minus_win(self.V0,i0,Tk,0.,pin) 
        wout = win + wowi                              # energy @ (f)

        popi = self.phiout_minus_phiin(fact,ga,r,i0,i1,Tk,0.,Tkp,0.,pin)
        pout = pin + popi                              # phase @ (f)
        
        dp = +(pout-self.PHOUT)    # delta phase @ (f)
        dw = +(wout-self.WOUT)     # delta energy @ (f)
        
        zf  = -dp*(c*be)/omeg
        zpf = ga/(ga+1)*dw/self.WOUT

        # f_particle = copy(self.particle)(tkin=wout)     # energy parameters from particle
        f_particle = copy(self.particle)(tkin=self.WOUT)  # energy parameters from soll
        gbf        = f_particle.gamma_beta

        fact = self.V0/(m0c2*gb*gbf)*i1
        gbi = gb
        if r > 0.:
            xp = gbi/gbf*xp-x/r*fact*Tk*sin(pin)
            yp = gbi/gbf*yp-y/r*fact*Tk*sin(pin)
        elif r == 0.:
            xp = gbi/gbf*xp
            yp = gbi/gbf*yp

        T = T + self.deltaW       # soll-energy gained
        f_track = NP.array([x,xp,y,yp,zf,zpf,T,1.,s,1.])

        if DEBUG_MAP == DEBUG_ON:    # for DEBUGGING
            dbTab1Row = [
                '{:8.4f}'.format(degrees(pout)),
                '{:8.4f}'.format(degrees(pin)),
                '{:8.4g}'.format(degrees(popi)),
                '{:8.4g}'.format(dp),
                '{:8.4g}'.format(wout),
                '{:8.4g}'.format(win),
                '{:8.4g}'.format(wowi),
                '{:8.4g}'.format(dw),
                '{:8.4g}'.format(self.WOUT),
                '{:8.3f}'.format(self.V0*1.e3),
                ]
            self.parent.dbTab1Rows.append(dbTab1Row)

            dbTab2Row = [
                '{:8.4g}'.format(x*1.e3),
                '{:8.4g}'.format(xp*1.e2),
                '{:8.4g}'.format(y*1.e3),
                '{:8.4g}'.format(yp*1.e3),
                '{:8.4g}'.format(z*1.e3),
                '{:8.4g}'.format(zp*1.e3),
                '{:8.4g}'.format(r*1.e3),
                '{:8.4g}'.format(Tk),
                '{:8.3g}'.format(Tkp),
                '{:8.4g}'.format(i0-1.),
                '{:8.4g}'.format(i1),
                ]
            self.parent.dbTab2Rows.append(dbTab2Row)
        return f_track
#todo: tests must be redone
def test0():
    import elements as ELM
    from tracks import Track
    
    print('-----------------------------------TEST 0----------------')
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    
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
    # track-point fields:          x   x'  y  y'  z   z'  T  1   s   1
    track = Track(start=NP.array([ x,  0., y, 0., z,  0., T, 1., 0., 1.]))
    ti = track.last()
    for i in range(1):
        DEBUG_TEST0('MAP:\n',track.last_str())
        tf = ttfg.map(ti)
        track.append(tf)
        DEBUG_TEST0('MAP:\n',track.last_str())
        ttfg.adjust_energy(tf[EKOO])    #enery adaptation
        ti = tf

def test1():
    import elements as ELM
    from tracks import Track
    
    print('-----------------------------------TEST 1----------------')
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    
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
    # start            x   x'  y  y'  z   z'  T  1   s   1
    start = NP.array([ x,  xp, y, yp, z,  zp, T, 1., 0., 1.])
    trck = Track(start=start)
    ti = trck.last()
    for i in range(anz+1):
        tf = ttfg.map(ti)
        trck.append(tf)
        z += delta
        start[4] = z
        ti = start
    DEBUG_TEST1('TRACK:\n',trck.asTable())

if __name__ == '__main__':
    test0()
    test1()