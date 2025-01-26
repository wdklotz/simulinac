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
from math import sin,cos,tan,sqrt,pi,degrees
from copy import copy
import numpy as NP

from setutil import FLAGS,PARAMS,Ktp,MDIM,Proton,DEBUG_ON,DEBUG_OFF
from Ez0 import SFdata
import elements as ELM

twopi = 2*pi
counter_of_polies = 0
trigger_poly_number = 12

class OXAL_G(ELM.RFG):
    """ OpenXAL RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf'], fieldtab=None):
        super().__init__(label, EzAvg, phisoll, gap, freq, SFdata=SFdata, particle=particle, position=position, aperture=aperture, dWf=dWf, mapping='oxal', fieldtab=fieldtab)
        if SFdata == None:
            raise RuntimeError('OXAL: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.map       = self.oxal_map   # OXAL's specific mapping method
            self.SFdata    = SFdata
            self.polies    = self.poly_slices(self.gap,self.SFdata)
            self.matrix    = self.make_matrix(self.polies,self.phisoll,self.particle)
            self.deltaW    = self.matrix[Ktp.T,Ktp.dT]
            self.particlef = Proton(particle.tkin + self.deltaW)

    def I0(self,k,d,cd,sd):
        return 2.*sd/k   # [cm]
    def I2(self,k,d,cd,sd):
        return 2.*(2.*d*cd/k**2+(d**2/k-2./k**3)*sd)   # [cm**3]
    def H1(self,k,d,cd,sd):
        return 2.*(sd/k**2-d*cd/k)     # [cm**2]
    def H3(self,k,d,cd,sd):
        return 2.*((3.*d**2/k**2-6./k**4)*sd-(d**3/k-6.*d/k**3)*cd) # [cm**4]
    def I4(self,k,d,cd,sd):
        return 2.*(d**4*sd/k-4./k*self.H3(k, d, cd, sd))   # [cm**5]
    def V0(self, poly):
        """ V0 A.Shishlo/J.Holmes (4.4.3) """
        # E0 = poly.E0                        # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)             # [cm]
        # v0 = v0*E0*self.dWf
        return v0                    # NOTE V0 is in [cm]
    def T(self, poly, k):
        """
        # poly = 1+a[1/cm]*z+b[1/cm**2]*z**2   
        # k = omega/(c*beta) [1/m]
        """
        """ T(k) A.Shishlo/J.Holmes (4.4.6) """
        b  = poly.b        # [1/cm**2]
        dz = poly.dz       # [cm]
        k  = k*1.e-2       # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        tk = (self.I0(k,dz,cd,sd) +b*self.I2(k,dz,cd,sd))  # [cm]
        tk = tk/v0         # []
        # if counter_of_polies == trigger_poly_number: DEBUG_ON('OXAL_slice:(T,k) {} [_]'.format((tk,k)))
        return tk   # []
    def S(self, poly, k):
        """ S(k) A.Shishlo/J.Holmes (4.4.7) """
        a  = poly.a
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        sk = a*self.H1(k,dz,cd,sd)    # [cm]
        sk = sk/v0                     # []
        # if counter_of_polies == trigger_poly_number: DEBUG_ON('OXAL_slice:(S,k) {} [_]'.format((sk,k)))
        return sk    # []
    def Tp(self, poly, k):
        """ 1st derivative T'(k) A.Shishlo/J.Holmes (4.4.8) """
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2      # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        tp = -(self.H1(k,dz,cd,sd) +b*self.H3(k,dz,cd,sd))  # [cm**2]
        tp = tp/v0   # [cm]
        # if counter_of_polies == trigger_poly_number: DEBUG_ON("OXAL_slice:(T',k) {} [cm]".format((tp,k)))
        return tp*1.e-2   #[m]
    def Sp(self, poly, k):
        """ 1st derivative S'(k) A.Shishlo/J.Holmes (4.4.9) """
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2      # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        sp = a* self.I2(k,dz,cd,sd)    #[cm**2]
        sp = sp/v0      # [cm]
        # if counter_of_polies == trigger_poly_number: DEBUG_ON("OXAL_slice:(S',k) {} [cm]".format((sp,k)))
        return sp*1.e-2      # [m]
    def Tpp(self, poly, k):
        """ 2nd derivative T''(k) """
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        tpp = -(self.I2(k,dz,cd,sd) + b*self.I4(k,dz,cd,sd) +a*self.H1(k,dz,cd,sd))  # [cm**3]
        tpp = tpp/v0    # [cm**2]
        # if counter_of_polies == trigger_poly_number: DEBUG_ON("OXAL_slice:(T'',k) {} [cm**2]".format((tpp,k)))
        return tpp*1.e-4      # [m**2]
    def Spp(self, poly, k):
        """ 2nd derivative S''(k) """
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        v0 = self.V0(poly) # [cm]
        cd = cos(k*dz)
        sd = sin(k*dz)
        spp = -a*self.H3(k,dz,cd,sd)   #[cm**3]
        spp = spp/v0   # [cm**2] 
        # if counter_of_polies == trigger_poly_number: DEBUG_ON("OXAL_slice:(S'',k) {} [cm**2]".format((spp,k)))
        return spp*1.e-4     # [m**2]
    def poly_slices(self, gap, SFdata):
        """Slice the RF gap"""
        slices = []
        zl = -gap/2.*100.   # [m] --> [cm]
        zr = -zl
        for poly in SFdata.EzPoly:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slices.append(poly)
        return slices
    def make_matrix(self, polies, phisoll, particle):
        c        = PARAMS['clight']
        m0c2     = particle.m0c2
        m0c3     = m0c2*c
        omega    = self.omega
        matrix   = NP.eye(MDIM,MDIM)

        # initialise loop variables
        p        = Proton(particle.tkin)
        phis     = phisoll
        for poly in polies:   # each poly is a slice of the full mat
            # global counter_of_polies   # debugging
            # counter_of_polies += 1     # debugging
            # z         = poly.dz*1.e-2     # [cm] ==> [m]
            Ws_in     = p.tkin
            betas_in  = p.beta
            gammas_in = p.gamma
            gbs_in    = p.gamma_beta
            gbs3_in   = gbs_in**3
            phis_in   = phis
            ks        = omega/(c*betas_in)

            # ptkin = p.tkin          # debugging
            # psdeg = degrees(phis)   # debugging
            
            qV0m   = self.V0(poly)*1.e-2*poly.E0     # NOTE change units V0 is in [cm]->[m]
            Tks    = self.T(poly,ks)
            Sks    = self.S(poly,ks)
            Tpks   = self.Tp(poly,ks)
            Spks   = self.Sp(poly,ks)
            Tppks  = self.Tpp(poly,ks)  
            Sppks  = self.Spp(poly,ks) 

            sphis  = sin(phis)
            cphis  = cos(phis)

            # kin.energy increase ref-particle
            DWs       = qV0m*(Tks*cphis - Sks*sphis)     # Shishlo 4.6.1
            Ws_out    = Ws_in + DWs
            # phase increase ref-particle
            Dphis = qV0m*omega/m0c3/gbs3_in*(Tpks*sphis + Spks*cphis)  # Shishlo 4.6.2
            phis_out   = phis_in + Dphis

            # OXAL-matrix 
            g3b2s_in    = gammas_in**3*betas_in**2    # gamma**3*beta**2 in
            gammas_out  = 1. + Ws_out/m0c2            # gamma  out
            gbs_out     = sqrt(gammas_out**2-1.)      # (gamma*beta) out
            betas_out   = gbs_out/gammas_out          # beta out
            g3b2s_out   = gammas_out**3*betas_out**2  # gamma-s**3*beta-s**2 out
            #======================================================================================="""        
            # (4.6.11) in Shishlo's paper:
            factor = qV0m*betas_out/m0c2/gbs3_in
            r44 = betas_out/betas_in + factor*omega/c/betas_in*(Tpks*cphis - Spks*sphis)
            r45 = factor*(3*gammas_in**2*(Tpks*sphis + Spks*cphis) + omega/c/betas_in*(Tppks*sphis+Sppks*cphis))
            # (4.6.9) in Shishlo's paper:
            r54 = qV0m*omega/(g3b2s_out*m0c3*betas_in)*(Tks*sphis + Sks*cphis)
            r55 = (g3b2s_in/g3b2s_out - qV0m*omega/(m0c3*betas_in*g3b2s_out)*(Tpks*cphis - Spks*sphis))
            #======================================================================================="""        
            # {z, Dbeta/betas}: linear sub matrix
            # NOTE: Shislo's Formeln sind fuer (z,Dbeta/betas) longitudinal
            mx = NP.eye(MDIM,MDIM)
            mx[Ktp.z, Ktp.z] = r44; mx[Ktp.z, Ktp.zp ] = r45
            mx[Ktp.zp,Ktp.z] = r54; mx[Ktp.zp, Ktp.zp] = r55
            # {x,x'}: linear sub-matrix
            factor = qV0m*omega/(2.*m0c3*gbs_out*gbs_in**2)
            mx[Ktp.xp,Ktp.x ] = -factor * (Tks*sphis + Sks*cphis)
            mx[Ktp.xp,Ktp.xp] = gbs_in/gbs_out
            # {y,y'}: linear sub-matrix
            mx[Ktp.yp,Ktp.y]  = mx[Ktp.xp, Ktp.x]
            mx[Ktp.yp,Ktp.yp] = mx[Ktp.xp, Ktp.xp]
            # energy and length increase
            mx[Ktp.T,Ktp.dT] = DWs
            mx[Ktp.S,Ktp.dS] = 0     # 0 length: oxal-gap is kick
            # add slice-matrix to oxal-matrix
            matrix = NP.dot(mx,matrix)

            # refresh loop variables
            # p = Proton(Ws_out)
            p = p(Ws_out)   # NOTE using the call method for Particle
            phis = phis_out
        return matrix
    def adjust_energy(self, tkin):
        adjusted = OXAL_G(self.label, self.EzAvg, self.phisoll, self.gap, self.freq, SFdata=self.SFdata, particle=Proton(tkin), position=self.position, aperture=self.aperture, dWf=self.dWf,fieldtab=self.fieldtab)
        return adjusted
    def oxal_map(self,i_track):
        gs2 = (self.particle.gamma)**2
        i_track = copy(i_track)
        i_track[Ktp.zp] = 1./gs2 * i_track[Ktp.zp]      # Dp/p -> Dbeta/beta
        f_track = NP.dot(self.matrix,i_track)
        gs2 = (self.particlef.gamma)**2
        f_track[Ktp.zp] = gs2 * f_track[Ktp.zp]         # Dbeta/beta -> Dp/p
        return f_track
if __name__ == '__main__':
    print('what?')
