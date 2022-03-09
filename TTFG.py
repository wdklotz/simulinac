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
import pprint, inspect

from setutil import PARAMS,I0,I1,tblprnt,arrprnt,Proton,FLAGS
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
from Ez0 import SFdata
import elements as ELM

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print(F'DEBUG_ON ==============>{file}: ', end="")
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

twopi          = 2*PI

class TTF_G(ELM.RFG):
    """Transition Time Factors RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247)"""
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(50.), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
        super().__init__(label, EzAvg, phisoll, gap, freq, particle, position, aperture, dWf, mapping='ttf')
        # TmStamp.stamp('OXAL init')
        if SFdata == None:
            raise RuntimeError('TTF_G: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.map       = self.ttf_g_map   # OXAL's specific mapping method
            self.SFdata    = SFdata
            self.polies    = self.poly_slices(self.gap,self.SFdata)
            # self.matrix    = self.make_matrix(self.polies,self.phisoll,self.particle)
            # self.deltaW    = self.matrix[Ktp.T,Ktp.dT]
            # self.particlef = Proton(particle.tkin + self.deltaW)

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
    def V0(self, poly):      # A.Shishlo/J.Holmes (4.4.3)
        # E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = (2*dz+2./3.*b*dz**3)             # [cm]
        # v0 = v0*E0*self.dWf
        return v0                             # NOTE [cm]
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
    def ttf_g_map(self, i_track):
        c          = PARAMS['clight']
        m0c2       = self.particle.m0c2
        m0c3       = m0c2*c
        omega      = self.omega

        # initialise loop variables
        p       = copy(self.particle)
        phis    = self.phisoll
        f_track = copy(i_track)
        for poly in self.polies:
            """ Map through this poly interval """
            x        = f_track[XKOO]       # [0]
            xp       = f_track[XPKOO]      # [1]
            y        = f_track[YKOO]       # [2]
            yp       = f_track[YPKOO]      # [3]
            z        = f_track[ZKOO]       # [4] z
            zp       = f_track[ZPKOO]      # [5] dp/p
            T        = f_track[EKOO]       # [6] kinetic energy ref
            S        = f_track[SKOO]       # [8] position
            # ref
            betas_in      = p.beta
            gammas_in     = p.gamma
            gbs_in        = p.gamma_beta
            gb3s_in       = gbs_in**3
            ks            = omega/(c*betas_in)
            Tk            = self.T(poly,ks)
            Tkp           = self.Tp(poly,ks)
            Sk            = self.S(poly,ks)
            Skp           = self.Sp(poly,ks)
            V0m           = self.V0(poly)*1.e-2        # NOTE V0 in [m]
            E0            = poly.E0
            phis_in       = phis                     
            Ws_in         = p.tkin
            cphis_in      = cos(phis_in)
            sphis_in      = sin(phis_in)
            """ Formel 4.3.1 A.Shishlo/J.Holmes """
            Ws_out_minus_Ws_in = V0m*E0*(Tk*cphis_in - Sk*sphis_in)
            Ws_out  = Ws_in + Ws_out_minus_Ws_in
            ps_out  = Proton(Ws_out)
            gammas_out = ps_out.gamma
            # tracked particle
            W_in          = zp*(gammas_in+1.)/gammas_in*Ws_in+Ws_in      # energy (i)  ~ (z')
            p_in          = Proton(W_in)
            gb_in         = p_in.gamma_beta
            r             = sqrt(x**2+y**2)            # radial coordinate
            K             = omega/(c*gbs_in)*r
            i0            = I0(K)                      # bessel function I0
            i1            = I1(K)                      # bessel function I1
            phi_in        = -z*omega/(c*betas_in)+phis_in # phase  (i)  ~ (-z)
            cphi_in       = cos(phi_in)
            sphi_in       = sin(phi_in)
            """ Formel 4.3.1 A.Shishlo/J.Holmes """
            W_out_minus_W_in = V0m*E0*i0*(Tk*cphi_in - Sk*sphi_in)
            W_out = Ws_in + W_out_minus_W_in
            p_out = Proton(W_out)
            gamma_out = p_out.gamma
            gb_out    = p.gamma_beta

            DW = W_out - Ws_out
            zp_out = gamma_out/(gamma_out+1.)*DW/Ws_out

            """ Formel 4.3.2 A.Shishlo/J.Holmes """
            faktor = V0m*E0*omega/m0c3/gb3s_in
            phis_out_minus_phis_in = faktor*(Tkp*cphis_in + Skp*sphis_in)
            phis_out = phis_in + phis_out_minus_phis_in

            """ Formel 4.3.2 A.Shishlo/J.Holmes """
            gamma_m = (gammas_out+gammas_in)/2.
            phi_out_minus_phi_in = faktor*(i0*(Tkp*cphi_in + Skp*sphi_in)+gamma_m*r*i1*(Tk*cphi_in+Sk*sphi_in))
            phi_out = phi_in + phi_out_minus_phi_in
                        
            z_out = gammas_out/(gammas_out+1)*DW/Ws_out       #TODO check w or DW/W
            """ Formel 4.3.3 A.Shishlo/J.Holmes """
            faktor = V0m*E0/(m0c2*gb_in*gb_out)*i1
            if r > 0.:
                xp = gb_in/gb_out*xp-x/r*faktor*(Tk*sphis_in + Sk*cphis_in)
                yp = gb_in/gb_out*yp-y/r*faktor*(Tk*sphis_in + Sk*cphis_in)
            elif r == 0.:
                xp = gb_in/gb_out*xp
                yp = gb_in/gb_out*yp

            # reset loop variables
            T = T + Ws_out_minus_Ws_in
            p    = ps_out
            phis = phis_out
            f_track = NP.array([x,xp,y,yp,z_out,zp_out,T,1.,S,1.])

        self.deltaW = Ws_out_minus_Ws_in
        self.particlef = ps_out
        return f_track
    def adjust_energy(self, tkin):
        adjusted = TTF_G(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted

if __name__ == '__main__':
    pass