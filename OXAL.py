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
import time

from setutil import PARAMS,DEBUG,DEBUG_ON,DEBUG_OFF,I0,I1,tblprnt,arrprnt,objprnt,Ktp,MDIM
from setutil import TmStamp
from Ez0 import SFdata

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
                # instantiate _TTF_Gslices
                slice = _OXAL_slice(parent, poly, particle)
                slices.append(slice)
            return slices

        #todo: need adjustment of phis to middle of gap ??
        def configure_slices(slices, phis, tkin):
            """adjust SOLL-energy dependence of slices"""
            TmStamp.stamp('config...')
            matrix = NP.eye(MDIM,MDIM)
            tkinIN = tkin
            phIN = phis
            for slice in slices:
                # call with time(aka phase) and energy @ entrance, return them @ exit
                tkin, phis = slice.adjust_slice_parameters(tkin,phis)
                # update Node matrix
                matrix = NP.dot(matrix,slice.matrix)
            deltaW   = tkin-tkinIN  # total energy advance
            deltaPhi = phis-phIN    # total phase advance
            return deltaW, deltaPhi, matrix

        TmStamp.stamp('OXAL init')
        # _OXAL attributes
        self.EzAvg    = parent.EzAvg
        self.gap      = parent.gap
        self.E0L      = self.EzAvg*self.gap
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.omega    = twopi*self.freq
        self.dWf      = parent.dWf
        self.SFdata   = parent.SFdata
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

            # configure slice for SOLL energy
            self._deltaW, dummy, matrix = configure_slices(self.slices, self.phis, self.tkin)

            # UPDATE linear NODE matrix
            parent.matrix = matrix

            # delayed  _OXAL attributes
            self._ttf = self._deltaW/(self.E0L*cos(self.phis)) if self.dWf == 1 else 1.
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
        # add SOLL-energy increase
        f_track[Ktp.T] += self._deltaW
        DEBUG_OFF('oxal-map ',f_track)
        return f_track

    def soll_map(self, i_track):
        TmStamp.stamp('soll_map')
        si,sm,sf = self.position
        f_track = copy(i_track)
        f_track[Ktp.T] += self._deltaW
        f_track[Ktp.S]  = sm
        DEBUG_OFF('oxal-soll ',f_track)
        return f_track
        
    def _full_gap_map(self, slices, i_track):
        """ The wrapper to slice mappings """
        TmStamp.stamp('full_gap_map')
        f_track = copy(i_track)
        for slice in slices:
            # map each slice
            f_track = slice.slice_map(f_track)
            pass
        return f_track
        
    def _T(self, poly, k):
        """ T(k) A.Shishlo/J.Holmes (4.4.6) """
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        # DEBUG_SLICE('_TTF_Gslice:_T: (T,k)',(t,k))
        return t

    def _S(self, poly, k):
        """ S(k) A.Shishlo/J.Holmes (4.4.7) """
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.-k*dz/tan(k*dz)
        s  = f1*f2
        # DEBUG_SLICE('_TTF_Gslice:_T: (T,k)',(t,k))
        return s

    def _Tp(self, poly, k):
        """ 1st derivative T'(k) A.Shishlo/J.Holmes (4.4.8) """
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp

    def _Sp(self, poly, k):
        """ 1st derivative S'(k) A.Shishlo/J.Holmes (4.4.9) """
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        sp  = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        sp  = sp*(dz**2-2./k**2+dz/tan(k*dz)*2/k)
        sp  = sp*1.e-2     # [cm] --> [m]
        return sp

    # def _Tpp(self, poly, k):
    #     """ 2nd derivative T''(k) """
    #     a   = poly.a
    #     b   = poly.b
    #     dz  = poly.dz
    #     k   = k*1.e-2      # [1/m] --> [1/cm]
    #     tpp = \
    #         2*(dz*k*(dz*k*(-6*b + k**2*(dz**2*b + 1))/tan(dz*k) + 6*b - k**2*(3*dz**2*b + 1))*cos(dz*k) 
    #         - (dz*k*(-6*b + k**2*(dz**2*b + 1))/tan(dz*k) + 6*b - k**2*(3*dz**2*b + 1))*sin(dz*k) 
    #         - (dz**2*k**2*(-6*b + k**2*(dz**2*b + 1))/sin(dz*k)**2 - 12*dz*b*k/tan(dz*k) + 18*b 
    #         - k**2*(3*dz**2*b + 1))*sin(dz*k))/(dz*k**5*(2./3.*dz**2*b + 2))
    #     return tpp

    #   def _Spp(self, poly, k):
    #     """ 2nd derivative S''(k) """
    #     a   = poly.a
    #     b   = poly.b
    #     dz  = poly.dz
    #     k   = k*1.e-2      # [1/m] --> [1/cm]
    #     spp = \
    #         2*a*(dz**3*k**3*cos(dz*k) - 3*dz**2*k**2*sin(dz*k)
    #         - 6*dz*k*cos(dz*k) + 6*sin(dz*k))/(dz*k**4*(2./3.*dz*b + 2))
    #     return spp

    def _V0(self, poly):
        """ V0 A.Shishlo/J.Holmes (4.4.3) """
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

# def dbg_slice(slice):
#     print('Win {:8.5} \t Phin {:8.3}'.format(slice.Wins,degrees(slice.phis)))

class _OXAL_slice(object):
    """ PyOrbit's openXAL linear RF-Gap Model """
    def __init__(self, parent, poly, particle):
        self.parent     = parent # the gap this slice belongs to
        self.particle   = copy(particle) # incoming SOLL particle
        # !!!ACHTUNG units!!! poly interval: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.poly       = poly 
        self.polydz     = poly.dz*1.e-2     # [cm] ==> [m]
        self.V0         = parent._V0(self.poly)
        self.phis       = None  # initialized in configure_slices
        self.Tks        = None  # initialized in adjust_slice_parameters
        self.Sks        = None  # initialized in adjust_slice_parameters
        self.Tpks       = None  # initialized in adjust_slice_parameters
        self.Spks       = None  # initialized in adjust_slice_parameters
        self.Tppks      = None  # initialized in adjust_slice_parameters   # not needed for linear model
        self.Sppks      = None  # initialized in adjust_slice_parameters   # not needed for linear model
        self.Wouts      = None  # initialized in adjust_slice_parameters
        self.dws        = None  # initialized in adjust_slice_parameters
        self.phis       = None  # initialized in adjust_slice_parameters
        self.particlef  = None  # initialized in adjust_slice_parameters

    def adjust_slice_parameters(self, tkin, phin):
        """ Adjust SOLL-energy dpendent parameters for this slice """
        self.particle(tkin)    # UPDATE tkin
        # TmStamp.stamp('adjust...')
        c      = PARAMS['lichtgeschwindigkeit']
        m0c2   = PARAMS['proton_mass']
        m0c3   = m0c2*c
        Wins   = self.particle.tkin
        betas  = self.particle.beta
        gammas = self.particle.gamma
        omega  = self.parent.omega
        phis   = phin
        qV0    = self.V0
        ks     = omega/(c*betas)
        
        Tks    = self.parent._T(self.poly,ks)
        Sks    = self.parent._S(self.poly,ks)
        Tpks   = self.parent._Tp(self.poly,ks)
        Spks   = self.parent._Sp(self.poly,ks)
        # Tppks  = self.parent._Tpp(self.poly,ks)   # not needed for linear model
        # Sppks  = self.parent._Spp(self.poly,ks)   # not needed for linear model
        sphis  = sin(phis)
        cphis  = cos(phis)

        # energy increase SOLL
        dws       = qV0*(Tks*cphis - Sks*sphis)
        Wouts     = Wins + dws
        particlef = copy(self.particle)(tkin=Wouts)
        # phase increase SOLL
        phiouts   = phis + omega*self.polydz/(c*betas)

        # oxal-matrix from SOLL aliases and variables
        gbs3_in     = 1./(gammas*betas)**3       # (gamma*beta)**3 in
        betas_in    = betas                      # beta in
        gammas_in   = gammas                     # gamma  in
        g3b2s_in    = gammas_in**3*betas_in**2   # gamma**3*beta**2 in
        gammas_out  = 1. + Wouts/m0c2            # gamma  out
        gbs_out     = sqrt(gammas_out**2-1)      # (gamma*beta) out
        betas_out   = gbs_out/gammas_out         # beta out
        g3b2s_out   = gammas_out**3*betas_out**2 # gamma-s**3*beta-s**2 out

#=======================================================================================
# (4.6.10) in Shishlo's paper mit SYMPY berechnet: 
# ACHTUNG! ist nicht dasselbe wie Shishlo's formel. Keine Tpp- & Spp-terme hier! 
        # DDPHI= # DDPHI= phi_out - phi_in 
        # -(Spks*sphis - Tpks*cphis)*dphi*gsbs3*omega*qV0/m0c3 
        # -(Spks*cphis + Tpks*sphis)*3*db2bs*omega*qV0/(betas**3*gammas*m0c3) 
        
        # - Sppks*cphis*db2bs**2*gsbs3*omega**2*qV0/(betas*c*m0c3)            O(2) (a) term mit gleichem factor in (4.6.10)
        # - Tppks*sphis*db2bs**2*gsbs3*omega**2*qV0/(betas*c*m0c3)            O(2) (b) term mit gleichem factor in (4.6.10)

        # - 3*Tpks*cphis*db2bs*dphi*omega*qV0/(betas**3*gammas*m0c3)          O(2)
        # + 3*Spks*db2bs*dphi*omega*qV0*sphis/(betas**3*gammas*m0c3)          O(2)
        # - Tppks*cphis*db2bs**2*dphi*gsbs3*omega**2*qV0/(betas*c*m0c3)       O(3)
        # + Sppks*db2bs**2*dphi*gsbs3*omega**2*qV0*sphis/(betas*c*m0c3)       O(3)
        # + 3*Sppks*cphis*db2bs**3*omega**2*qV0/(betas**4*c*gammas*m0c3)      O(3)
        # + 3*Tppks*db2bs**3*omega**2*qV0*sphis/(betas**4*c*gammas*m0c3)      O(3)
        # - 3*Sppks*db2bs**3*dphi*omega**2*qV0*sphis/(betas**4*c*gammas*m0c3) O(4)
        # + 3*Tppks*cphis*db2bs**3*dphi*omega**2*qV0/(betas**4*c*gammas*m0c3) O(4)
# Facit: Formeln (4.6.10) und (4.6.11) im paper sind nicht korreckt. Terme (a) & (b) sind O(2) und nicht O(1) !
# Facit: die sympy-Rechnung ist genauer als die Formeln im paper !!
#=======================================================================================        
# (4.6.9) in Shishlo's paper:
        # dbeta2beta_out = db2bs*(g3b2s_in/g3b2s_out-qV0*omega/(m0c3*betas_in*g3b2s_out)*(Tpks*cphis-Spks*sphis)) 
        # + z*qV0*omega/(g3b2s_out*m0c3*betas_in)*(Tks*sphis+Sks*cphis)
# (4.6.11) in Shishlo's paper:
        # z_out = betas_out/betas_in*z 
        # + qV0*betas_out/(m0c2*gbs3_in)*((3*gammas_in**2*(Tpks*sphis+Spks*cphis)
        # + omega/(c*betas_in)*(Tppks*sphis+Sppks*cphis))*db2bs 
        # + omega/(c*betas_in)*(Tpks*cphis-Spks*sphis)*z)
#=======================================================================================        

        # {z, delta-beta/betas}: linear sub matrix
        mx = NP.eye(MDIM,MDIM)
        # implementation (4.6.9) und (4.6.11) from Shishlo's paper
        mx[Ktp.z,Ktp.z ] = betas_out/betas_in + qV0*betas_out/(m0c2*gbs3_in)*omega/(c*betas_in)*(Tpks*cphis-Spks*sphis)
        # mx[Ktp.z,Ktp.zp] = qV0*betas_out/(m0c2*gbs3_in)*(3*gammas_in**2*(Tpks*sphis+Spks*cphis)+omega/(c*betas_in)*(Tppks*sphis+Sppks*cphis))
        # (4.6.11) mit meiner sympy Korrektur: keine Tpp- und Spp-terme in (4.6.11)
        mx[Ktp.z,Ktp.zp] = qV0*betas_out/(m0c2*gbs3_in)*(3*gammas_in**2*(Tpks*sphis+Spks*cphis))
        mx[Ktp.zp,Ktp.z] = qV0*omega/(g3b2s_out*m0c3*betas_in)*(Tks*sphis+Sks*cphis)
        mx[Ktp.zp,Ktp.zp]= (g3b2s_in/g3b2s_out - qV0*omega/(m0c3*betas_in*g3b2s_out)*(Tpks*cphis-Spks*sphis))

        # {x,x',y,y'}: linear sub-matrix
        gbs_in     = gammas*betas
        gbs_out    = sqrt(gammas_out**2-1.)
        facxy      = -qV0*omega/(2.*m0c3*gbs_out*gbs_in**2)
        mx[Ktp.xp,Ktp.x]  = facxy * (Tks*sphis + Sks*cphis) # mx(x',x) * xin
        mx[Ktp.xp,Ktp.xp] = gbs_in/gbs_out                  # mx(x',x')* xin'
        mx[Ktp.yp,Ktp.y]  = mx[Ktp.xp,Ktp.x]                # mx(y',y) * yin
        mx[Ktp.yp,Ktp.yp] = mx[Ktp.xp,Ktp.xp]               # mx(y',y')* yin'

        # energy and length increase
        mx[Ktp.T,Ktp.dT] = dws
        mx[Ktp.S,Ktp.dS] = 0     # oxal-gap is kick of DKD - no length increase
        
        # _OXAL_slice attributes needed by slice_map
        self.matrix     = mx
        self.gammas     = gammas
        self.gammas_out = gammas_out
        return Wouts, phiouts

    def slice_map(self, i_track):
        """Map through this slice from position (i) to (f)"""
        # TmStamp.stamp('slice_map')
        # z      = i_track[Ktp.z]       # [4] z~(phi-phis)
        zp       = i_track[Ktp.zp]      # [5] delta-p/p
        track    = copy(i_track)
        # local aliases
        gammas     = self.gammas
        gammas_out = self.gammas_out        

        db2bs = DbetaToBetaFromDp2p(gammas,zp)      # delta-beta/betas in
        track[Ktp.zp]   = db2bs                     # zp ==> db2bs
        f_track = NP.dot(self.matrix,track)         # matrix
        zp_out  = gammas_out**2*f_track[Ktp.zp]     # db2bs ==> zp
        f_track[Ktp.zp] = zp_out
        DEBUG_OFF(repr(self.__dict__))
        DEBUG_OFF('oxal-slice ',f_track)
        # dbg_slice(self)
        return f_track

def test0():
    from bunch import Tpoint, Track
    from elements import RFC,RFG
    
    print('-----------------------------------TEST 0----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']*1.8055 # [Mv/m] EzPeak/EzAvg fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,EzPeak)
    
    oxal = RFG(gap=0.048,SFdata=SF_tab,mapping='oxal')
    tkin = 50.
    oxal.adjust_energy(tkin=tkin)
    if DEBUG_TEST0():
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
    from bunch import Tpoint, Track
    from elements import RFC,RFG
    
    print('-----------------------------------TEST 1----------------')
    input_file='SF_WDK2g44.TBL'
    EzPeak = PARAMS['EzAvg']
    SF_tab = SFdata(input_file,EzPeak)
    
    oxal = RFG(gap=0.048,SFdata=SF_tab,mapping='oxal')
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
