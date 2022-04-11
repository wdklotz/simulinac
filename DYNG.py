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
import numpy as NP
from copy import copy
from math import degrees,radians,sqrt,cos,sin,pi,pow
# from functools import partial
# import warnings
from collections import namedtuple
import unittest

import elements as ELM
# from setutil import DEB, arrprnt, PARAMS, tblprnt, Ktp, WConverter
from setutil import DEB, PARAMS, Proton, FLAGS
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Polyval

DEBUG_OFF = DEB.get('OFF')
DEBUG_ON  = DEB.get('ON')

EzAvg_test =1.    # 1 Mev constant gap field

class DYN_G(ELM.RFG): 
    """ DYNAC's RF-gap model. Numerical computations in an accelerating gap or in a cavity; E.TANKE and S.VALERO 3-Oct-2016 """
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
        super().__init__(label, EzAvg, phisoll, gap, freq, particle, position, aperture, dWf, mapping='dyn')
        if SFdata == None:
            raise RuntimeError('DYNG: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.SFdata    = SFdata
            self.map       = self.dynac_map   # OXAL's specific mapping method
            self.polies    = self.poly_slices(self.gap,self.SFdata)
            self.particlef = None
            self.ttf       = 0.

    def Ez0t(self,z,t,omega,phi,EzAvg=EzAvg_test):
        return EzAvg*cos(omega*t+phi)
    def dEz0tdt(self,z,t,omega,phi,EzAvg=EzAvg_test):
        return -omega*EzAvg*sin(omega*t+phi)
    def Panofsky_test(self,gap,beta,lamb,phis,EzAvg=EzAvg_test):
        """ calculate energy increase with Panofsky's formula """
        x = gap/beta/lamb
        ttf = NP.sinc(x)
        E0TL = EzAvg*ttf*gap
        DW = E0TL*cos(phis)
        return DW
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
    def Picht(self,g,DgDz,x,xp,y,yp,back=False):
            """ A Picht transformation from cartesian to reduced variables """
            g2m1     = g**2-1.
            # g2m1p05  = sqrt(g2m1)
            g2m1p025 = pow(g2m1,0.25)
            if not back:
                # r=(x,y) rp=(xp,yp)
                X  = g2m1p025*x
                Y  = g2m1p025*y
                XP = g2m1p025*xp+0.5*x*g*DgDz/g2m1
                YP = g2m1p025*yp+0.5*y*g*DgDz/g2m1
                return (X,XP,Y,YP)
            else:
                # r=(X,Y) rp=(Xp,Yp)
                X  = x; XP = xp; Y = y; YP = yp
                x  = X/g2m1p025
                y  = Y/g2m1p025
                xp = (XP-0.5*X*g*DgDz/g2m1)/g2m1p025
                yp = (YP-0.5*Y*g*DgDz/g2m1)/g2m1p025
                return (x,xp,y,yp)
    def dynac_map(self,i_track): 
        #TODO: use non_const beta, i.e. gamma[i], i=1,4,1 in inegrals
        #TODO: find position for multiplication with dWf-flag
        #TODO: use equivalent E-field instead of SFdata
        """ Mapping from (i) to (f) """  
        # stpfac = self.stpfac
        c      = PARAMS['clight']
        m0c2   = PARAMS['proton_mass']
        m0c3   = m0c2*c
        omega  = self.omega
        lamb   = self.lamb
        """ initialize loop variables """
        p           = copy(self.particle)
        Ws_in       = p.tkin
        #TODO self.deltaW = -Ws_in
        phis     = self.phisoll    # phase @ gap middle
        f_track  = copy(i_track)
        for poly in self.polies:    
            """ Map through a single poly interval 
            For each poly interval a numerical step-by-step method based on the 5 point Boole's will be applied. The step size in the azimuthally direction z is divided in 4 parts equivalent 
            lengths."""     
            # Azimutally positions z0,z1,z2,z3,z4; step size h
            h  = (poly.zr-poly.zl)*1.e-2     # stepsize [m]
            h4 = h/4.
            z0 = poly.zl*1.e-2     # [m]
            z1 = z0 + h4
            z2 = z1 + h4
            z3 = z2 + h4
            z4 = z3 + h4
            
            g0        = p.gamma             # ref @ in
            b0        = p.beta              # ref @ in
            W0        = p.tkin
            gb0       = g0*b0
            gb03      = gb0**3
            K1        = omega**2/4/c**2/gb03
            b0c       = b0*c
            h4b0c     = h/4/b0c
            t0        = z0/b0c
            t1        = t0 + h4b0c
            t2        = t1 + h4b0c
            t3        = t2 + h4b0c
            t4        = t3 + h4b0c

            # Ez0   = self.SFdata.Ez0t(z0*100.,t0,omega,phis)
            # Ez1   = self.SFdata.Ez0t(z1*100.,t1,omega,phis)
            # Ez2   = self.SFdata.Ez0t(z2*100.,t2,omega,phis)
            # Ez3   = self.SFdata.Ez0t(z3*100.,t3,omega,phis)
            # Ez4   = self.SFdata.Ez0t(z4*100.,t4,omega,phis)
            Ez0   = self.Ez0t(z0*100.,t0,omega,phis)
            Ez1   = self.Ez0t(z1*100.,t1,omega,phis)
            Ez2   = self.Ez0t(z2*100.,t2,omega,phis)
            Ez3   = self.Ez0t(z3*100.,t3,omega,phis)
            Ez4   = self.Ez0t(z4*100.,t4,omega,phis)
            # phi0  = degrees(omega*t0+phis)      # phase @ pos z0
            # phi1  = degrees(omega*t1+phis)
            # phi2  = degrees(omega*t2+phis)
            # phi3  = degrees(omega*t3+phis)
            # phi4  = degrees(omega*t4+phis)
            # phisdegs = degrees(phis)
            I1    = h/90.*(7.*Ez0 + 32.*Ez1 + 12.*Ez2 + 32.*Ez3 + 7.*Ez4)
            I2    = h**2/90.*(8.*Ez1 + 6.*Ez2 + 24.* Ez3 + 7.*Ez4)

            """ Energy """
            # DWs = self.Panofsky_test(h, b0, lamb, omega*t2+phis)
            # Dgs_pan   = DWs/m0c2
            Dgs       = I1/m0c2
            DgDz      = Dgs/h

            """ Picht """
            x         = f_track[XKOO]       # [0]
            xp        = f_track[XPKOO]      # [1]
            y         = f_track[YKOO]       # [2]
            yp        = f_track[YPKOO]      # [3]
            xr,xpr,yr,ypr \
                      = self.Picht(g0,DgDz,x,xp,y,yp,back=False)
            R02       = xr**2 + yr**2
            R0R0p     = xr*xpr + yr*ypr

            """ Energy """
            g4s  = g0 + Dgs        # s @ pos z4
            Dgp  =  ((1+R02*K1)*I1 + R0R0p*K1*I2)/m0c2
            g4p  = g0 + Dgp        # p @ z4

            zp_  = f_track[ZPKOO]   # dp/p @ z0
            Dg0p = (g0+1)/g0 * zp_   # p @ z0
            Dg4p = Dg0p+Dgp-Dgs     # p @ z4
            zp   = g0/(g0+1) * Dg4p # dp/p @ z4 TODO g0 or g4s?
            
            Ws_out = (g4s-1) * m0c2
            # Wp_out = (g4p-1) * m0c2
            ps_out = Proton(Ws_out)
            # p_out  = Proton(Wp_out)
            # gbs3   = ps_out.gamma_beta**3
            # gbp3   = p_out.gamma_beta**3

            """ Phase """
            # I3 = h**2/90./gb03*(8.*Ez1 + 6.*Ez2 + 24.*Ez3 + 7.*Ez4)
            # I4 = h**3/90./gb03*(2.*Ez1 + 3.*Ez2 + 18.*Ez3 + 7.*Ez4)
            # I3_ = h**2/90.*(8.*Ez1 + 6.*Ez2 + 24.*Ez3 + 7.*Ez4)
            # I4_ = h**3/90.*(2.*Ez1 + 3.*Ez2 + 18.*Ez3 + 7.*Ez4)
            # I3s = I3_/gbs3
            # I3p = I3_/gbp3
            # I4p = I4_/gbp3

            # Dts  = 1/m0c3*I3s
            # Dts_  = h/b0c
            # Dtp  = 1/m0c3*((1+R02*K1)*I3p + R0R0p*K1*I4p)
            # t4s  = t0 + Dts + h/b0c
            # t4p  = t0 + Dtp + h/b0c
            # bs0 = b0
            # bs4 = ps_out.beta
            z_  = f_track[ZKOO]   
            # z__ = - h/b0c * z_
            # z   = bs4/bs0 * z__
            z = z_ + b0*h/gb03*(g0-1)/g0 * zp

            # Dphi = (- z_/b0c + Dtp - Dts) * omega
            # W0p = W0*(1+(g0+1)/g0*zp_)
            # b0p = Proton(W0p).beta
            # t4s = t0 + h/b0c + Dts
            # t4p = t0 - z_/b0p/c + h/b0p/c + Dtp
            # Dphi = (t4p-t4s) * omega
            # z    = - Dphi*b0c/omega

            """ transverse x,xp,y,yp """
            g2m1 = g0**2-1
            G1 = 1./2/m0c3*pow(g2m1,-1.5)
            # Ez0p = self.SFdata.dEz0tdt(z0*100,t0,omega)
            # Ez1p = self.SFdata.dEz0tdt(z1*100,t1,omega)
            # Ez2p = self.SFdata.dEz0tdt(z2*100,t2,omega)
            # Ez3p = self.SFdata.dEz0tdt(z3*100,t3,omega)
            # Ez4p = self.SFdata.dEz0tdt(z4*100,t4,omega)
            Ez0p = self.dEz0tdt(z0*100,t0,omega,phis)
            Ez1p = self.dEz0tdt(z1*100,t1,omega,phis)
            Ez2p = self.dEz0tdt(z2*100,t2,omega,phis)
            Ez3p = self.dEz0tdt(z3*100,t3,omega,phis)
            Ez4p = self.dEz0tdt(z4*100,t4,omega,phis)
            J1 = h/90.*G1*(7.*Ez0p + 32.*Ez1p + 12.*Ez2p + 32.*Ez3p + 7.*Ez4p)
            J2 = h**2/90.*G1*(8.*Ez1p + 6.*Ez2p + 24.*Ez3p + 7.*Ez4p)

            Dxpr     = xr*J1 + xpr*J2
            xpr_out  = xpr + Dxpr
            Dypr     = yr*J1 + ypr*J2
            ypr_out  = ypr + Dypr

            J3 = h**3/90.*G1*((2.*Ez1p + 3.*Ez2p + 18.*Ez3p + 7.*Ez4p))

            Dxr     = xr*J2 + xpr*J3
            xr_out  = xr + Dxr + h*xpr
            Dyr     = yr*J2 + ypr*J3
            yr_out  = yr + Dyr + h*ypr

            x,xp,y,yp = self.Picht(g4s,DgDz,xr_out,xpr_out,yr_out,ypr_out,back=True)

            """ Reset der loop Variablen """
            p        = ps_out
            T        = f_track[EKOO]
            S        = f_track[SKOO]
            f_track  = NP.array([x,xp,y,yp,z,zp,T,1,S,1])

        # leaving loop poly intervals
        self.particlef = p
        self.deltaW    = Ws_out - Ws_in
        self.ttf       = self.deltaW/EzAvg_test/h/cos(phis)/len(self.polies)
        f_track[EKOO]  = Ws_out
        DEBUG_OFF('dyn-map {}'.format(i_track))
        DEBUG_OFF('dyn-map {}'.format(f_track))
        return f_track

    def adjust_energy(self, tkin):
        adjusted = DYN_G(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted
class TestDynacGapMapping(unittest.TestCase):
    def test_DYNG_REF_mapping(self):
        print('----------------------------------test_DYNG_REF_mapping')
        label    = 'dyn_gap_test'
        phisoll  = radians(-25.)
        gap      = 0.022
        freq     = 800e6
        tkin     = 50
        particle = Proton(tkin)
        dWf      = 1
        fname    = 'Superfish/SF_WDK2g44.TBL'
        gap_cm   = gap*100     # Watch out!
        EzPeak   = 2.2
        Ezdata   = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
        EzAvg    = Ezdata.EzAvg
        dyngap   = DYN_G(label,EzAvg,phisoll,gap,freq,SFdata=Ezdata,particle=particle,dWf=dWf)

        i_track  = NP.array([0.,0.,0.,0.,0.,0.,tkin,1,0,1])
        i_track  = NP.array([1e-3, 0., 1e-3, 0.,0,0.,tkin,1,0,1])
        i_track  = NP.array([0., 0., 0., 0.,0,+1e-4,tkin,1,0,1])
        f_track  = dyngap.map(i_track)
if __name__ == '__main__':
    unittest.main()
