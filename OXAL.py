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
import sys
import numpy as NP
import unittest
import elements as ELM

from math import sin,cos,tan,sqrt,pi,degrees,radians
from copy import copy
from setutil import FLAGS,PARAMS,Ktp,MDIM,Proton,DEBUG_ON,DEBUG_OFF,Proton
from Ez0 import SFdata

twopi = 2*pi
counter_of_polies = 0
trigger_poly_number = 12

class OXAL_G(ELM.RFG):
    """ OpenXAL RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self, label, EzPeak, phisoll, cavlen, freq, sfdata, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
        super().__init__(label, EzPeak, phisoll, 0., cavlen,freq, sfdata, particle, position, aperture, dWf)
        # self.viseo     = 0.25
        # self.label     = label
        # self.EzPeak    = EzPeak
        # self.phisoll   = phisoll
        self.cavlen    = cavlen
        self.length    = cavlen
        # self.freq      = freq
        # self.SFdata    = sfdata
        # self.particle  = particle
        # self.position  = position
        # self.aperture  = aperture
        # self.dWf       = dWf
        self.omega     = twopi*self.freq
        self.polies    = self.poly_slices()
        self.matrix    = self.make_matrix()
        self.deltaW    = self.matrix[Ktp.T,Ktp.dT]
        self.particlef = Proton(particle.tkin + self.deltaW)

    def V0(self, poly):
        """ V0 A.Shishlo/J.Holmes (4.4.3) """
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = E0*(2*dz+2./3.*b*dz**3)*1.e-2    # [MV]
        return v0 
    def T(self, poly, k):    # A.Shishlo/J.Holmes (4.4.6)
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        DEBUG_OFF('TTF_G: (T,k) {}'.format((t,k)))
        return t
    def S(self, poly, k):    # A.Shishlo/J.Holmes (4.4.7)
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.-k*dz/tan(k*dz)
        s  = f1*f2
        DEBUG_OFF('TTF_G: (S,k) {}'.format((s,k)))
        return s
    def Tp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp
    def Sp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.9)
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
        return 0
    def Spp(self, poly, k):
        """ 2nd derivative S''(k) """
        return 0
    def poly_slices(self):
        """Slice the RF cavity"""
        L = self.cavlen
        sfdata = self.SFdata
        slices = []
        zl = -L/2.*100.   # [m] --> [cm]
        zr = -zl
        for poly in sfdata.polies:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slices.append(poly)
        return slices
    def make_matrix(self):
        polies   = self.polies
        c        = PARAMS['clight']
        m0c2     = self.particle.m0c2
        m0c3     = m0c2*c
        omega    = self.omega
        matrix   = NP.eye(MDIM,MDIM)

        # initialise loop variables
        psoll        = Proton(self.particle.tkin)
        phis     = self.phisoll

        for poly in polies:   # each poly is a slice of the full mat
            # global counter_of_polies   # debugging
            # counter_of_polies += 1     # debugging
            # z = poly.dz*1.e-2     # [cm] ==> [m]

            # IN variables
            Ws_in     = psoll.tkin
            betas_in  = psoll.beta
            gammas_in = psoll.gamma
            gbs_in    = psoll.gamma_beta
            gbs3_in   = gbs_in**3
            phis_in   = phis
            g3b2s_in  = gammas_in**3*betas_in**2    # gamma**3*beta**2 in
            g2s_in    = gammas_in**2
            ks        = omega/(c*betas_in)

            # ptkin = psoll.tkin      # debugging
            # psdeg = degrees(phis)   # debugging
            
            qV0    = self.V0(poly)     # [MV]
            Tks    = self.T(poly,ks)
            Sks    = self.S(poly,ks)
            Tpks   = self.Tp(poly,ks)
            Spks   = self.Sp(poly,ks)
            Tppks  = self.Tpp(poly,ks)  
            Sppks  = self.Spp(poly,ks) 

            sphis  = sin(phis)
            cphis  = cos(phis)

            # kin.energy increase ref-particle
            DWs       = qV0*(Tks*cphis - Sks*sphis)     # Shishlo 4.6.1
            # phase increase ref-particle
            Dphis = qV0*omega/m0c3/gbs3_in*(Tpks*sphis + Spks*cphis)  # Shishlo 4.6.2
            Ws_out    = Ws_in + DWs
            phis_out  = phis_in + Dphis

            # OUT variables
            gammas_out  = 1. + Ws_out/m0c2            # gamma  out
            gbs_out     = sqrt(gammas_out**2-1.)      # (gamma*beta) out
            betas_out   = gbs_out/gammas_out          # beta out
            g3b2s_out   = gammas_out**3*betas_out**2  # gamma-s**3*beta-s**2 out
            g2s_out     = gammas_out**2
            #======================================================================================="""        
            # OXAL-matrix 
            # (4.6.11) in Shishlo's paper:
            factor  = qV0*betas_out/m0c2/gbs3_in
            factor1 = qV0*omega/m0c3/betas_in/g3b2s_in
            r44 = betas_out/betas_in + factor*omega/c/betas_in*(Tpks*cphis - Spks*sphis)
            r45 = factor*(3*gammas_in**2*(Tpks*sphis + Spks*cphis) + omega/c/betas_in*(Tppks*sphis+Sppks*cphis))
            # (4.6.9) in Shishlo's paper:
            r54 = factor1*(Tks*sphis + Sks*cphis)
            r55 = (g3b2s_in/g3b2s_out - factor1*(Tpks*cphis - Spks*sphis))
            #======================================================================================="""        
            # {z, Dbeta/betas}: linear sub matrix
            # NOTE: Shislo's Formeln sind fuer (z,Dbeta/betas) longitudinal
            mx = NP.eye(MDIM,MDIM)
            mx[Ktp.z, Ktp.z] = r44;         mx[Ktp.z, Ktp.zp ] = r45/g2s_in  # apply conversion DBeta/Beta=gamma**(-2)*Dp/p
            mx[Ktp.zp,Ktp.z] = r54*g2s_out; mx[Ktp.zp, Ktp.zp] = r55         # apply conversion Dp/p=gamma**2*DBeta/Beta
            # {x,x'}: linear sub-matrix
            factor2 = qV0*omega/(2.*m0c3*gbs_out*gbs_in**2)
            mx[Ktp.xp,Ktp.x ] = -factor2 * (Tks*sphis + Sks*cphis)
            mx[Ktp.xp,Ktp.xp] = gbs_in/gbs_out
            # {y,y'}: linear sub-matrix
            mx[Ktp.yp,Ktp.y]  = mx[Ktp.xp, Ktp.x]
            mx[Ktp.yp,Ktp.yp] = mx[Ktp.xp, Ktp.xp]
            # energy and length increase
            mx[Ktp.T,Ktp.dT] = DWs
            mx[Ktp.S,Ktp.dS] = 0     # 0 length: oxal-gap is kick
            # left muöltiplication of slice-matrix to oxal-matrix
            matrix = NP.dot(mx,matrix)

            # refresh loop variables
            psoll = psoll(Ws_out)   # NOTE using the call method for Particle
            phis = phis_out
        return matrix
    def adjust_energy(self, tkin):
        adjusted = OXAL_G(self.label, self.EzPeak, self.phisoll, self.cavlen, self.freq, self.SFdata, particle=Proton(tkin), position=self.position, aperture=self.aperture, dWf=self.dWf)
        return adjusted

class TestOxalEnergyMapping(unittest.TestCase):
    def test_OXAL(self):
        """ testing the OXAL mapping for acceleration """
        injection_energies = [6,12,24,50,100,150,200]
        EzPeak = 2.0; phisoll = radians(-30.); gap = 0.044; freq = 750.e6; fieldtab="SF/SF_WDK2g44.TBL"
        gap_cm = 100.*gap   # gap in cm
        sfdata = SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,IntgInterval=gap_cm)
        EzAvg = sfdata.EzAvg
        ID = "OXAL"
        for injection_energy in injection_energies:
            oxal = OXAL_G(ID,EzAvg,phisoll,gap,freq,SFdata=sfdata,particle=Proton(injection_energy))
            tvectori = NP.array([0, 0, 0, 0, 0, 0, injection_energy, 1, 0, 1])
            tvectoro = NP.dot(oxal.matrix,tvectori)
            print(f'EzPeak={EzPeak}, gap={gap}, freq={freq*1e-6}, Wi={tvectori[6]:3.3f}, Wo={tvectoro[6]:3.3f}, dW={oxal.matrix[6,7]:.3e}')
        return

if __name__ == '__main__':
    unittest.main()
