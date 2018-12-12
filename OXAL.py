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
from math import sin,cos,tan,radians,degrees,sqrt,pi
from copy import copy
import numpy as NP

from setutil import PARAMS,DEBUG,I0,I1,tblprnt,arrprnt,objprnt,Ktp
import elements as ELM
from Ez0 import SFdata

# DEBUG__*
def DEBUG_ON(string,arg='',end='\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg='',end='\n'):
    pass
    
DEBUG_SLICE = DEBUG_OFF
DEBUG_TEST0 = DEBUG_ON
DEBUG_TEST1 = DEBUG_ON

twopi          = 2*pi

class _OXAL(object):
    """ OpenXAL RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247)"""
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
                slice = _OXAL_slice(parent, poly, particle)
                slices.append(slice)
            return slices

        #todo: need adjustment of phis to middle of gap ??
        def configure_slices(slices, phis, tkin):
            """adjust energy dependence of slices"""
            tkinIN = tkin
            phIN = phis
            for slice in slices:
                # call with time(aka phase) and energy @ entrance, return them @ exit
                tkin, phis = slice.adjust_slice_parameters(tkin,phis)
            deltaW   = tkin-tkinIN  # total energy advance
            deltaPhi = phis-phIN    # total phase advance
            return deltaW, deltaPhi

        # _OXAL
        self.EzAvg    = parent.EzAvg
        self.gap      = parent.gap
        self.E0L      = self.EzAvg*self.gap
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.omega    = twopi*self.freq
        self.dWf      = parent.dWf
        self.SFdata   = parent.SFdata
        self.matrix   = parent.matrix
        self.particle = parent.particle
        self.tkin     = self.particle.tkin
        self.position = parent.position
        self._deltaW    = None # initailized by configure_slices()
        self._particlef = None # initailized by configure_slices()
        if parent.SFdata == None:
            raise RuntimeError('_TTF_G: missing E(z) table - STOP')
            sys.exit(1)
        else:
             # slice the gap
            self.slices = make_slices(self, self.gap, self.SFdata, self.particle)
            # slice energy dependence
            self._deltaW, dummy = configure_slices(self.slices, self.phis, self.tkin)
            self._ttf = self._deltaW/(self.E0L*cos(self.phis)) if self.dWf == 1 else 1.
            # UPDATE linear NODE matrix with deltaW
            self.matrix[Ktp.T,Ktp.dT] = self._deltaW
            self._particlef = copy(self.particle)(self.particle.tkin + self._deltaW)
            parent['slices'] = self.slices   # satify test0()
            # for slice in self.slices:
            #     print(repr(slice.__dict__))
                # dbg_slice(slice)
            pass

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
        # full map through sliced openXAL gap-model
        f_track = self._full_gap_map(self.slices, f_track)
        f_track[Ktp.T] += self._deltaW
        DEBUG_OFF('oxal-map ',f_track)
        return f_track

    def soll_map(self, i_track):
        si,sm,sf = self.position
        f_track = copy(i_track)
        f_track[Ktp.T] += self._deltaW
        f_track[Ktp.S]  = sm
        DEBUG_OFF('oxal-soll ',f_track)
        return f_track
        
    def _full_gap_map(self, slices, i_track):
        """ The wrapper to slice mappings """
        f_track = copy(i_track)
        for slice in slices:
            # map each slice with openXAL gap-model
            f_track = slice.slice_map(f_track)
            # relativistic scaling. Is it needed?
            # z = f_track[ZKOO]
            # betai = self.particle.beta
            # tkin  = self.particle.tkin
            # betaf = Proton(tkin=tkin+self.deltaW).beta
            # z = betaf/betai*z
            # f_track[ZKOO] = z
            pass
        return f_track
        
    def _T(self, poly, k):    # A.Shishlo/J.Holmes (4.4.6)
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        # DEBUG_SLICE('_TTF_Gslice:_T: (T,k)',(t,k))
        return t

    def _S(self, poly, k):    # A.Shishlo/J.Holmes (4.4.7)
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.-k*dz/tan(k*dz)
        s  = f1*f2
        # DEBUG_SLICE('_TTF_Gslice:_T: (T,k)',(t,k))
        return s

    def _Tp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp

    def _Sp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.9)
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        sp  = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        sp  = sp*(dz**2-2./k**2+dz/tan(k*dz)*2/k)
        sp  = sp*1.e-2     # [cm] --> [m]
        return sp

    def _Tpp(self, poly, k):
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tpp = \
            2*(dz*k*(dz*k*(-6*b + k**2*(dz**2*b + 1))/tan(dz*k) + 6*b - k**2*(3*dz**2*b + 1))*cos(dz*k) 
            - (dz*k*(-6*b + k**2*(dz**2*b + 1))/tan(dz*k) + 6*b - k**2*(3*dz**2*b + 1))*sin(dz*k) 
            - (dz**2*k**2*(-6*b + k**2*(dz**2*b + 1))/sin(dz*k)**2 - 12*dz*b*k/tan(dz*k) + 18*b 
            - k**2*(3*dz**2*b + 1))*sin(dz*k))/(dz*k**5*(2./3.*dz**2*b + 2))
        return tpp

    def _Spp(self, poly, k):
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        spp = \
            2*a*(dz**3*k**3*cos(dz*k) - 3*dz**2*k**2*sin(dz*k)
            - 6*dz*k*cos(dz*k) + 6*sin(dz*k))/(dz*k**4*(2./3.*dz*b + 2))
        return spp

    def _V0(self, poly):    # A.Shishlo/J.Holmes (4.4.3)
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)*1.e-2       # [cm] --> [m]
        v0 = v0*E0*self.dWf
        return v0

def DbetaToBetaFromDp2p(gamma,Dp2p):
    return Dp2p/gamma**2
def DphiFromZ(omega,c,beta,z):
    return -omega/(c*beta)*z
def ZFromDphi(omega,c,beta,dphi):
    return -(c*beta)/omega*dphi
def Dw2wFromDp2p(gamma,Dp2p):
    return (gamma+1)/gamma*Dp2p
def Dp2pFromDw2w(gamma,Dw2w):
    return gamma/(gamma+1)*Dw2w
def dbg_slice(slice):
    print('Win {:8.5} \t Phin {:8.3}'.format(slice.Wins,degrees(slice.phis)))

class _OXAL_slice(object):
    """ PyOrbit's openXAL RF-Gap Model """
    def __init__(self, parent, poly, particle):
        self.parent     = parent # the gap this slice belongs to
        self.particle   = copy(particle) # incoming SOLL particle
        # !!!ACHTUNG units!!! poly interval: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.poly       = poly 
        self.length     = poly.dz*1.e-2     # [cm] ==> [m]
        self.V0         = parent._V0(self.poly)
        self.phis       = None  # initialized in configure_slices
        # self.ks         = None  # initialized in adjust_slice_parameters
        self.Tks        = None  # initialized in adjust_slice_parameters
        self.Sks        = None  # initialized in adjust_slice_parameters
        self.Tpks       = None  # initialized in adjust_slice_parameters
        self.Spks       = None  # initialized in adjust_slice_parameters
        self.Tppks      = None  # initialized in adjust_slice_parameters
        self.Sppks      = None  # initialized in adjust_slice_parameters
        self.Wouts      = None  # initialized in adjust_slice_parameters
        self.dws        = None  # initialized in adjust_slice_parameters
        self.phis       = None  # initialized in adjust_slice_parameters
        self.particlef  = None  # initialized in adjust_slice_parameters

    def adjust_slice_parameters(self, tkin, phi):
        """ Adjust energy-dpendent SOLL parameters for this slice """
        self.particle(tkin)    # UPDATE tkin
        c      = PARAMS['lichtgeschwindigkeit']
        m0c2   = PARAMS['proton_mass']
        m0c3   = m0c2*c
        Wins   = self.particle.tkin
        betas  = self.particle.beta
        gammas = self.particle.gamma
        omega  = self.parent.omega
        phis   = phi
        qV0    = self.V0
        ks     = omega/(c*betas)
        
        Tks    = self.parent._T(self.poly,ks)
        Sks    = self.parent._S(self.poly,ks)
        Tpks   = self.parent._Tp(self.poly,ks)
        Spks   = self.parent._Sp(self.poly,ks)
        # Tppks  = self.parent._Tpp(self.poly,ks)
        # Sppks  = self.parent._Spp(self.poly,ks)
        sphis  = sin(phis)
        cphis  = cos(phis)

        # energy
        dws       = qV0*(Tks*cphis - Sks*sphis)
        Wouts     = Wins + dws
        particlef = copy(self.particle)(tkin=Wouts)
        # phase
        Phiouts   = phis + omega*self.length/(c*betas)

        self.mx = NP.eye(ELM.MDIM,ELM.MDIM)
        # tranverse: linear submatrix matrix {x,x',y,y'}
        gbINs    = gammas*betas
        gammaOUTs= 1 + Wouts/m0c2
        gbOUTs   = sqrt(gammaOUTs**2-1.)
        facxy    = -qV0*omega/(2.*m0c3*gbOUTs*gbINs**2)
        self.mx[Ktp.xp,Ktp.x]  = facxy * (Tks*sphis + Sks*cphis) # mx(x',x) * xin
        self.mx[Ktp.xp,Ktp.xp] = gbINs/gbOUTs                    # mx(x',x')* xin'
        self.mx[Ktp.yp,Ktp.y]  = facxy * (Tks*sphis + Sks*cphis) # mx(y',y) * yin
        self.mx[Ktp.yp,Ktp.yp] = gbINs/gbOUTs                    # mx(y',y')* yin'

        # energy and length increase
        self.mx[Ktp.T,Ktp.dT] = dws
        self.mx[Ktp.S,Ktp.dS] = self.length
        
        # variables needed by slice_map
        self.Wins  = Wins
        self.phis  = phis
        self.betas = betas
        self.gammas= gammas
        self.omega = omega
        self.qV0   = qV0
        self.Tks   = Tks
        self.Sks   = Sks
        self.Tpks  = Tpks
        self.Spks  = Spks
        # self.Tppks = Tppks
        # self.Sppks = Sppks
        self.cphis = cphis
        self.sphis = sphis
        self.dws   = dws
        self.Wouts = Wouts
        self.gsbs3= 1./(gammas*betas)**3
        return Wouts, Phiouts
#todo: check again when to use (bg)-IN and (bg)-OUT 
#todo: check against Shishlo's formulas
    def slice_map(self, i_track):
        """Map through this slice from position (i) to (f)"""
        z        = i_track[Ktp.z]       # [4] z~(phi-phis)
        zp       = i_track[Ktp.zp]      # [5] delta-p/p
        track = copy(i_track)
        # locals
        c      = PARAMS['lichtgeschwindigkeit']
        m0c2   = PARAMS['proton_mass']
        m0c3   = m0c2*c
        betas  = self.betas
        gammas = self.gammas
        gsbs3 = self.gsbs3
        omega  = self.omega
        qV0    = self.qV0
        Wins   = self.Wins
        # phis   = self.phis
        Tks    = self.Tks
        Sks    = self.Sks
        Tpks   = self.Tpks
        Spks   = self.Spks
        # Tppks  = self.Tppks
        # Sppks  = self.Sppks
        cphis  = self.cphis
        sphis  = self.sphis
        dws    = self.dws                                       # kin. energy increase in gap SOLL
        Wouts  = self.Wouts
        
        Win   = Wins + Dw2wFromDp2p(gammas,zp)*Wins             # W @ IN
        db2bs = DbetaToBetaFromDp2p(gammas,zp)                  # delta-beta/betas
        dphi  = DphiFromZ(omega,c,betas,z)                      # delta-phi

        # phase increase SOLL        
        DPHIS = gsbs3*omega*qV0*(Spks*cphis + Tpks*sphis)/m0c3

        # phase increase PARTICLE
        DPHI = \
        + gsbs3*omega*qV0/m0c3 *                       (Spks*cphis + Tpks*sphis)
        - dphi*gsbs3*omega*qV0/m0c3 *                  (Spks*sphis + Tpks*cphis)
        - 3 * db2bs*omega*qV0/(betas**3*gammas*m0c3) * (Spks*cphis - Tpks*sphis)
        # formula made with symy. O(2) and higher terms commented
        # dphi + phis 
        # + 3*Tppks*db2bs**3*omega**2*qV0*sphis/(betas**4*c*gammas*m0c3) 
        # + 3*Spks*db2bs*dphi*omega*qV0*sphis/(betas**3*gammas*m0c3) 
        # - Sppks*cphis*db2bs**2*gsbs3*omega**2*qV0/(betas*c*m0c3) 
        # + Sppks*db2bs**2*dphi*gsbs3*omega**2*qV0*sphis/(betas*c*m0c3) 
        # + 3*Sppks*cphis*db2bs**3*omega**2*qV0/(betas**4*c*gammas*m0c3) 
        # - 3*Sppks*db2bs**3*dphi*omega**2*qV0*sphis/(betas**4*c*gammas*m0c3) 
        # - 3*Tpks*cphis*db2bs*dphi*omega*qV0/(betas**3*gammas*m0c3) 
        # - Tppks*cphis*db2bs**2*dphi*gsbs3*omega**2*qV0/(betas*c*m0c3) 
        # - Tppks*db2bs**2*gsbs3*omega**2*qV0*sphis/(betas*c*m0c3) 
        # + 3*Tppks*cphis*db2bs**3*dphi*omega**2*qV0/(betas**4*c*gammas*m0c3) 
        
        # kin. energy increase in gap PARTICLE
        dw = \
        qV0 * dphi *                (-Sks*cphis - Tks*sphis) 
        qV0 *                       (-Sks*sphis + Tks*cphis)
        db2bs*omega*qV0/(betas*c) * (+Spks*sphis - Tpks*cphis)
        # formula made with symy. O(2) and higher terms commented
        # + Spks*cphis*db2bs*dphi*omega*qV0/(betas*c) 
        # + Tpks*db2bs*dphi*omega*qV0*sphis/(betas*c) 
        
        Wouts = Wins + dws           # SOLL
        Wout  = Win + dw             # PARTICLE
        DWout = Wout - Wouts         # delta-delta
        gammaouts = 1.+ Wouts/m0c2
        gbouts    = sqrt(gammaouts**2-1.)
        betaouts  = gbouts/gammaouts
        Dp2pout   = Dp2pFromDw2w(gammaouts,DWout/Wouts)
        
        dphio     = dphi+DPHI-DPHIS
        zpout     = Dp2pout
        zout      = ZFromDphi(omega,c,betaouts,dphio)
        
        f_track = NP.dot(self.mx,track)
        f_track[Ktp.z]  = zout
        f_track[Ktp.zp] = zpout
        DEBUG_OFF(repr(self.__dict__))
        DEBUG_OFF('oxal-slice ',f_track)
        # dbg_slice(self)
        return f_track

def test0():
    import elements as ELM
    from bunch import Tpoint, Track
    
    print('-----------------------------------TEST 0----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']*1.8055 # [Mv/m] EzPeak/EzAvg fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,EzPeak)
    
    oxal = ELM.RFG(gap=0.048,SFdata=SF_tab,mapping='oxal')
    tkin = 50.
    oxal.adjust_energy(tkin=tkin)
    if DEBUG_TEST0 == DEBUG_ON:
        print('TTFG: oxal.__dict__',oxal.__dict__)      # for DEBUGGING
        slices = oxal['slices']
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
        tf = oxal.map(ti())
        tpf = Tpoint(tf)
        track.addpoint(tpf)
        DEBUG_TEST0('MAP:\n',track.getpoints()[-1].as_str())
        oxal.adjust_energy(tf[Ktp.T])    #enery adaptation
        ti = tpf

def test1():
    import elements as ELM
    from bunch import Tpoint, Track
    
    print('-----------------------------------TEST 1----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']
    SF_tab = SFdata(input_file,EzPeak)
    
    oxal = ELM.RFG(gap=0.048,SFdata=SF_tab,mapping='oxal')
    tkin = 150.
    oxal.adjust_energy(tkin=tkin)
    DEBUG_TEST1('TTFG: oxal.__dict__',oxal.__dict__)      # for DEBUGGING
    
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
        tf = oxal.map(ti())
        tpf = Tpoint(tf)
        track.addpoint(tpf)
        z += delta
        tf[4] = z
        ti = Tpoint(tf)
    DEBUG_TEST1('TRACK-POINTS:\n',track.as_table())

if __name__ == '__main__':
    test0()
    test1()