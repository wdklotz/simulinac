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
from math import sin,cos,tan,sqrt,pi
from copy import copy
import numpy as NP
# import time
import pprint, inspect

from setutil import PARAMS,Ktp,MDIM
from Ez0 import SFdata
import elements as ELM

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

twopi = 2*pi

class OXAL(ELM.RFG):
    """ OpenXAL RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(50.), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
        super().__init__(label, EzAvg, phisoll, gap, freq, SFdata, particle, position, aperture, dWf, mapping='oxal')
        # TmStamp.stamp('OXAL init')
        if self.SFdata == None:
            raise RuntimeError('OXAL: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.matrix = NP.eye(MDIM,MDIM)
            polies      = self.poly_slices(self.gap,self.SFdata)
            self.make_mx(polies,self.phisoll,self.particle)
            self.deltaW = self.matrix[Ktp.T,Ktp.dT]
            self.particlef = Proton(particle.tkin + self._deltaW)

    def V0(self, poly):
        """ V0 A.Shishlo/J.Holmes (4.4.3) """
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)*1.e-2       # [cm] --> [m]
        v0 = v0*E0*self.dWf
        return v0
    def T(self, poly, k):
        """ T(k) A.Shishlo/J.Holmes (4.4.6) """
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        DEBUG_SLICE('_TTF_Gslice:_T: (T,k) {}'.format((t,k)))
        return t
    def S(self, poly, k):
        """ S(k) A.Shishlo/J.Holmes (4.4.7) """
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.-k*dz/tan(k*dz)
        s  = f1*f2
        DEBUG_SLICE('_TTF_Gslice:_T: (T,k) {}'.format((s,k)))
        return s
    def Tp(self, poly, k):
        """ 1st derivative T'(k) A.Shishlo/J.Holmes (4.4.8) """
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp
    def Sp(self, poly, k):
        """ 1st derivative S'(k) A.Shishlo/J.Holmes (4.4.9) """
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        sp  = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        sp  = sp*(dz**2-2./k**2+dz/tan(k*dz)*2/k)
        sp  = sp*1.e-2     # [cm] --> [m]
        return sp
    def Tpp(self, poly, k):
        """ 2nd derivative T''(k) """
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
    def Spp(self, poly, k):
        """ 2nd derivative S''(k) """
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        spp = \
            2*a*(dz**3*k**3*cos(dz*k) - 3*dz**2*k**2*sin(dz*k)
            - 6*dz*k*cos(dz*k) + 6*sin(dz*k))/(dz*k**4*(2./3.*dz*b + 2))
        return spp    
    def poly_slices(self, gap, SFdata):
        """Slice the RF gap"""
        slices = []
        zl = -gap/2.*100.   # [m] --> [cm]
        zr = -zl
        #    E0z = 0.
        #    z = 0.
        for poly in SFdata.EzPoly:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            # instantiate _TTF_Gslices
            # slice = OXAL_slice(parent, poly, particle)
            slices.append(poly)
        return slices
    def make_mx(self, polies, phisoll, particle):
        c        = PARAMS['clight']
        m0c2     = Proton.m0c2
        m0c3     = m0c2*c
        omega    = self.omega

        # initialise loop variables
        p        = particle
        phis     = phisoll
        for poly in polies:   # each poly is a slice of the full matrix
            polydz = poly.dz*1.e-2     # [cm] ==> [m]
            Wins   = p.tkin
            betas  = p.beta
            gammas = p.gamma
            ks     = omega/(c*betas)
            
            qV0    = self.V0(poly)
            Tks    = self.T(self.poly,ks)
            Sks    = self.S(self.poly,ks)
            Tpks   = self.Tp(self.poly,ks)
            Spks   = self.Sp(self.poly,ks)
            # Tppks  = self.parent._Tpp(self.poly,ks)   # not needed for linear model
            # Sppks  = self.parent._Spp(self.poly,ks)   # not needed for linear model

            sphis  = sin(phis)
            cphis  = cos(phis)

            # energy increase ref-particle
            dws       = qV0*(Tks*cphis - Sks*sphis)
            Wouts     = Wins + dws
            # phase increase SOLL
            phiouts   = phis + omega*polydz/(c*betas)

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

            # add slice-matrix to oxal-matrix
            self.matrix = NP.dot(mx,self.matrix)

            # refresh loop variables
            p = Proton(Wouts)
            phis = phiouts
    def adjust_energy(self, tkin):
        adjusted = OXAL(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted
