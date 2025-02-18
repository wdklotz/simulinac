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
import Igap
from math import sin,cos,tan,radians,degrees,sqrt
# from copy import copy
import numpy as NP
import unittest
from setutil import PARAMS,I0,I1,tblprnt,arrprnt,Proton,FLAGS
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,DSKOO
from setutil import DEBUG_ON,DEBUG_OFF
from Ez0 import SFdata
# import elements as ELM

twopi = 2.*pi

class TTF_G(Igap.Igap):
    """Transition Time Factors RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247)"""
    # def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf'],fieldtab=None):
    #     super().__init__(label, EzAvg, phisoll, gap, freq, SFdata=SFdata, particle=particle, position=position, aperture=aperture, dWf=dWf, mapping='ttf', fieldtab=fieldtab)
    def __init__(self):
        pass

    def configure(self,**kwargs):
        self.length       = 0. # because it's a kick
        self.dWf          = FLAGS['dWf']
        self.mapping      = 'ttf'        # map model
        self.kwargs       = kwargs
        self.label        = 'TTF' 

        self.EzPeak    = kwargs.get('EzPeak',None)
        self.phisoll   = kwargs.get('phisoll',None)
        self.cavlen    = kwargs.get('cavlen',None)
        self.freq      = kwargs.get('freq',None)
        self.SFdata    = kwargs.get('SFdata',None)
        self.particle  = kwargs.get('particle',Proton(50.))
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)

        self.lamb      = PARAMS['clight']/self.freq
        self.gap       = self.cavlen*0.57   # 57%: a best guess ?
        self.omega     = twopi*self.freq
        self.polies    = self.poly_slices(gap,SFdata)
        self.ttf       = 0
        self.matrix    = None
        self.deltaW    = None
        self.particlef = None
        self.master    = None

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return self.ttf_map(i_track)

    def toString(self):
        return 'base mapping: TTF_M.ttf_map'

    def isAccelerating(self):
        return True

    def waccept(self): pass

    def register_mapper(self,master):
        master.register_mapping(self)
        pass

    def accept_register(self,master):
        self.master = master
        pass

    def adjust_energy(self,tkin): pass

    def T3D_matrix(self):
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m       = NP.eye(MDIM,MDIM)
        E0L     = self.EzPeak*self.gap
        qE0LT   = E0L*self.ttf
        deltaW  = E0L*self.ttf*cos(self.phisoll)
        Wavg    = self.particle.tkin+self.deltaW/2.   # average tkin
        pavg    = Proton(Wavg)
        bavg    = pavg.beta
        gavg    = pavg.gamma
        m0c2    = pavg.e0
        kz      = twopi*E0L*self.ttf*sin(self.phisoll)/(m0c2*bavg*bavg*self.lamb)
        ky      = kx = -0.5*kz/(gavg*gavg)
        bgi     = self.particle.gamma_beta
        bgf     = self.particlef.gamma_beta
        bgi2bgf = bgi/bgf
        m       = NP.eye(MDIM,MDIM)
        m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
        m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
        m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
        m[EKOO, DEKOO] = deltaW
        m[SKOO, DSKOO]  = 0.
        self.matrix = m
        return


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
        """ Slice the RF gap """
        slices = []
        zl = -gap/2.*100.   # [m] --> [cm]
        zr = -zl
        for poly in SFdata.EzPoly:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slices.append(poly)
        return slices
    
    def ttf_map(self, i_track):
        c          = PARAMS['clight']
        m0c2       = self.particle.m0c2
        m0c3       = m0c2*c
        omega      = self.omega

        """ initialise loop variables """
        p       = self.particle
        phis    = self.phisoll
        self.deltaW = 0.
        for poly in self.polies:
            """ Map through this poly interval.
            Das Sollteilchen Formel 4.3.1 und 4.3.2 A.Shishlo/J.Holmes """
            gammas_in     = p.gamma
            betas_in      = p.beta
            Ws_in         = p.tkin 
            ks            = omega/(c*betas_in)
            Tk            = self.T(poly,ks)
            self.ttf      = self.ttf + Tk    # Shishlo's ttf of each poly  (4.4.4)
            Sk            = self.S(poly,ks)
            L0            = self.V0(poly)*1.e-2  # NOTE V0 in [m]
            qE0L          = poly.E0*L0
            phis_in       = phis                     
            cphis_in      = cos(phis_in)
            sphis_in      = sin(phis_in)
            Ws_out_minus_Ws_in = qE0L*(Tk*cphis_in - Sk*sphis_in)    # 4.3.1

            """ Sollenergie Out """
            Ws_out        = Ws_in + Ws_out_minus_Ws_in
            ps_out        = Proton(Ws_out)

            gbs_in        = p.gamma_beta
            gb3s_in       = gbs_in**3
            Tkp           = self.Tp(poly,ks)
            Skp           = self.Sp(poly,ks)
            faktor        = qE0L*omega/m0c3/gb3s_in
            phis_out_minus_phis_in = faktor*(Tkp*cphis_in + Skp*sphis_in)   # 4.3.2

            """ Sollphase Out """
            phis_out      = phis_in + phis_out_minus_phis_in

            """ Offteilchen Koordinaten IN """
            x         = i_track[XKOO]       # [0]
            xp        = i_track[XPKOO]      # [1] Dx/Ds
            y         = i_track[YKOO]       # [2]
            yp        = i_track[YPKOO]      # [3] Dy/Ds
            z         = i_track[ZKOO]       # [4] z
            zp        = i_track[ZPKOO]      # [5] dp/p

            max_r     = 0.05                # max radial excursion [m]
            r         = sqrt(x**2+y**2)     # radial coordinate
            if r > max_r:
                raise UTIL.OutOfRadialBoundEx(S)
            Kr        = omega/(c*gbs_in)*r
            i0        = I0(Kr)                      # bessel function I0 
            i1        = I1(Kr)                      # bessel function I1

            """ Delta-W = (gamma+1)/gamma * Delta-p/p * W """
            DW_in     = (gammas_in+1)/gammas_in * zp  * Ws_in 
            """ Delta-phi = -360/(beta*lambda) * z """
            Dphi_in   = -z * omega/(betas_in*c) # die z-Koordinate des Offteilchens als Dphi [rad]  
            phi_in    = Dphi_in + phis_in 
            cphi_in   = cos(phi_in)
            sphi_in   = sin(phi_in)

            """ Offteilchen Energiedifferenz """
            W_out_minus_W_in = qE0L*i0*(Tk*cphi_in - Sk*sphi_in)    # 4.3.1
            W_out     = Ws_in + W_out_minus_W_in
            p_out     = Proton(W_out)

            """ Mittelwert gamma_m """
            gammas_out = ps_out.gamma
            gamma_m    = (gammas_out+gammas_in)/2.

            """ Offteilchen Phasendifferenz """
            phi_out_minus_phi_in = faktor*(i0*(Tkp*cphi_in + Skp*sphi_in)+gamma_m*r*i1*(Tk*cphi_in+Sk*sphi_in))   # 4.3.2
            phi_out    = phi_in + phi_out_minus_phi_in

            """ transversale Koordinaten OUT Polyintervalls (4.3.3 Shishlo/Holmes) """
            gbs_out = ps_out.gamma_beta
            faktor  = qE0L/(m0c2*gbs_in*gbs_out)*i1
            if r > 0.:
                xp = gbs_in/gbs_out*xp-x/r*faktor*(Tk*sphis_in + Sk*cphis_in)
                yp = gbs_in/gbs_out*yp-y/r*faktor*(Tk*sphis_in + Sk*cphis_in)
            elif r == 0.:
                xp = gbs_in/gbs_out*xp
                yp = gbs_in/gbs_out*yp

            """ longitudinalen Koordinaten OUT Polyintervalls (4.3.3 Shishlo/Holmes) """
            betas_out = ps_out.beta
            """ z: Delta-phi = -360/(beta*lambda) * z """
            z_out     = - betas_out*c/omega * (Dphi_in + phi_out_minus_phi_in - phis_out_minus_phis_in)
            """ Delta-p/p: Delta-W = (gamma+1)/gamma * Delta-p/p * W """
            zp_out    = gammas_out/(gammas_out+1) * (DW_in + W_out_minus_W_in - Ws_out_minus_Ws_in)/Ws_out
   
            """ Offteilchens OUT """
            T = f_track[EKOO]       # [6] kinetic energy ref
            S = f_track[SKOO]       # [8] position
            T = Ws_out
            f_track = NP.array([x,xp,y,yp,z_out,zp_out,T,1.,S,1.])

            """ Reset der loop Variablen """
            p    = ps_out
            phis = phis_out

        self.particlef = ps_out
        self.deltaW    = ps_out.tkin
        self.ttf       = self.ttf/len(self.polies) # gap's ttf  (better as Panofski?)
        return f_track

class TestTransitTimeFactorsGapModel(unittest.TestCase):
    def test_TTFG_mapping(self):
        print('----------------------------------test_TTFG_mapping')
        phisoll  = radians(-25.)
        gap      = 0.040
        freq     = 800e6
        particle = Proton(50.)
        dWf      = 1
        fname    = 'SF/SF_WDK2g44.TBL'
        gap_cm   = gap*100     # Watch out!
        EzPeak   = 10.0
        Ezdata   = SFdata.field_data(fname,EzPeak=EzPeak,gap=gap_cm)
        EzAvg    = Ezdata.EzAvg
        ttfg     = TTF_G("rfg-test",EzAvg,phisoll,gap,freq,SFdata=Ezdata,particle=particle,dWf=dWf)

        i_track  = NP.array([0,0,0,0,0,0,50,1,0,1])
        f_track  = ttfg.map(i_track)
        # print(i_track); print(f_track)
        for i in range(len(f_track)):
            self.assertAlmostEqual(f_track[i],NP.array([0,0,0,0,0,0,50.13019,1,0,1])[i],msg="f_track",delta=1e-4)

        i_track  = NP.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 50, 1, 0, 1])
        f_track  = ttfg.map(i_track)
        # print(i_track); print(f_track)
        for i in range(len(f_track)):
            self.assertAlmostEqual(f_track[i],NP.array([1e-3, 1.0136e-3, 1e-3, 1.0136e-3, 1.0011e-3, 0.96408e-3, 50.13019, 1, 0, 1])[i],msg="f_track",delta=1e-4)
if __name__ == '__main__':
    unittest.main()