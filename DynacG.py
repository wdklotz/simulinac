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
import math as MATH
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

class DYNG(ELM.RFG):
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
    def Picht(self,gamma,x,xp,y,yp,back=False):
            """ A Picht transformation from cartesian to reduced variables """
            g2m1           = gamma**2-1.
            sqrt_g2m1      = MATH.sqrt(g2m1)
            sqrt_sqrt_g2m1 = MATH.sqrt(sqrt_g2m1)
            if not back:
                # r=(x,y) rp=(xp,yp)
                X = sqrt_sqrt_g2m1*x
                Y = sqrt_sqrt_g2m1*y
                XP = sqrt_sqrt_g2m1*xp+0.5*X*gamma/g2m1
                YP = sqrt_sqrt_g2m1*yp+0.5*y*gamma/g2m1
                return (X,XP,Y,YP)
            else:
                # r=(X,Y) rp=(Xp,Yp)
                X = x; XP=xp; Y=y; YP=yp
                x = X/sqrt_sqrt_g2m1
                y = Y/sqrt_sqrt_g2m1
                xp = (XP-0.5*gamma*X/g2m1)/sqrt_sqrt_g2m1
                yp = (YP-0.5*gamma*Y/g2m1)/sqrt_sqrt_g2m1
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
        """ initialize loop variables """
        p           = copy(self.particle)
        Ws_in       = p.tkin
        self.deltaW = -Ws_in
        phis_in     = self.phisoll
        f_track     = copy(i_track)
        for poly in self.polies:    
            """ Map through a single poly interval 
            For each poly interval a numerical step-by-step method based on the 5 point Boole's will be applied. The step size in the azimuthally direction z is divided in 4 parts equivalent 
            lengths."""     
            # Azimutally positions z0,z1,z2,z3,z4; step size h
            h  = (poly.zr-poly.zl)*1.e-2     # stepsize [m]
            z0 = poly.zl*1.e-2     # [m]
            z1 = z0 + h/4.
            z2 = z0 + h/2.
            z3 = z0 + 3.*h/4.
            z4 = z0 + h
            
            gammas_in = p.gamma             # ref @ in
            gbs_in    = p.gamma_beta
            gbs3_in   = gbs_in**3
            K1        = omega**2/4/c**2/gbs3_in
            x         = f_track[XKOO]       # [0]
            xp        = f_track[XPKOO]      # [1]
            y         = f_track[YKOO]       # [2]
            yp        = f_track[YPKOO]      # [3]
            Rx_in,Rxp_in,Ry_in,Ryp_in = self.Picht(gammas_in,x,xp,y,yp,back=False)
            betas_in  = p.beta
            bc        = betas_in*c
            t0        = z0/bc
            t1        = t0 + h/(4*bc)
            t2        = t0 + h/(2*bc)
            t3        = t0 + 3.*h/(4*bc)
            t4        = t0 + h/bc

            Ez0   = self.SFdata.Ez0t(z0*100.,t0,omega,phis_in)
            Ez1   = self.SFdata.Ez0t(z1*100.,t1,omega,phis_in)
            Ez2   = self.SFdata.Ez0t(z2*100.,t2,omega,phis_in)
            Ez3   = self.SFdata.Ez0t(z3*100.,t3,omega,phis_in)
            Ez4   = self.SFdata.Ez0t(z4*100.,t4,omega,phis_in)
            I1    = h/90.*(7.*Ez0 + 32.*Ez1 + 12.*Ez2 + 32.*Ez3 + 7.*Ez4)
            I2    = h**2/90.*(8.*Ez1 + 6.*Ez2 + 24.* Ez3 + 7.*Ez4)
            R02   = Rx_in**2 + Ry_in**2
            R0R0p = Rx_in*Rxp_in + Ry_in*Ryp_in

            """ Energy """
            Dgammas_out = I1/m0c2
            gammas_out  = gammas_in + Dgammas_out        # ref @ out
            zp          = f_track[ZPKOO]                 # dp/p
            gamma_in    = (gammas_in+1)/gammas_in * zp   # dp/p -> DW/W
            Dgamma      = 1./m0c2*((1+R02*K1)*I1 + R0R0p*K1*I2)
            gamma_out   = gamma_in + Dgamma

            I3 = h**2/90./gbs3_in*(8.*Ez1 + 6.*Ez2 + 24.*Ez3 + 7.*Ez4)
            I4 = h**3/90./gbs3_in*(2.*Ez1 + 3.*Ez2 + 18.*Ez3 + 7.*Ez4)

            """ Phase """
            Dts         = I3/m0c2
            t4s_out     = t0 + Dts + h/bc           # time ref out
            Dphis_out   = omega * t4s_out
            phis_out    = phis_in + Dphis_out       # ref @ out    
            Dt          = 1/m0c3*((1+R02*K1)*I3 + R0R0p*K1*I4)
            t4_out      = t0 + Dt + h/bc
            Dphi_out    = omega * t4_out           
            z           = f_track[ZKOO]             # z
            phi_in      = -z * omega/bc             # z -> Dphi/phi
            phi_out     = phi_in + Dphi_out         # phase @ out

            """ transverse x,xp,y,yp """
            g2m1 = gammas_in-1
            G1 = 1./2/m0c3*MATH.sqrt(1./g2m1**3)
            Ez0p = self.SFdata.dEz0tdt(z0*100,t0,omega,phis_in)
            Ez1p = self.SFdata.dEz0tdt(z1*100,t1,omega,phis_in)
            Ez2p = self.SFdata.dEz0tdt(z2*100,t2,omega,phis_in)
            Ez3p = self.SFdata.dEz0tdt(z3*100,t3,omega,phis_in)
            Ez4p = self.SFdata.dEz0tdt(z4*100,t4,omega,phis_in)
            J1 = h/90.*G1*(7.*Ez0p + 32.*Ez1p + 12.*Ez2p + 32.*Ez3p + 7.*Ez4p)
            J2 = h**2/90.*G1*(8.*Ez1p + 6.*Ez2p + 24.*Ez3p + 7.*Ez4p)

            DRxp_out = Rx_in*J1 + Rxp_in*J2
            Rxp_out  = Rxp_in + DRxp_out
            DRyp_out = Ry_in*J1 + Ryp_in*J2
            Ryp_out  = Ryp_in + DRyp_out

            J3 = h**3/90./G1*((2.*Ez1p + 3.*Ez2p + 18.*Ez3p + 7.*Ez4p))

            DRx_out = Rx_in*J2 + Rxp_in*J3
            Rx_out = Rx_in + DRx_out + h*Rxp_in
            DRy_out = Ry_in*J2 + Ryp_in*J3
            Ry_out = Ry_in + DRy_out + h*Ryp_in

            x,xp,y,yp = self.Picht(gammas_out,Rx_out,Rxp_out,Ry_out,Ryp_out,back=True)

            """ Reset der loop Variablen """
            p       = Proton(gammas_out*m0c2)
            phis_in = phis_out
            T       = f_track[EKOO]
            S       = f_track[SKOO]
            f_track     = NP.array([x,xp,y,yp,z,zp,T,1,S,1])

        # leaving loop poly intervals
        self.particlef = p
        self.deltaW    = self.deltaW + p.tkin
        self.ttf       = self.ttf/len(self.polies)
        f_track[EKOO]  = Ws_in + self.deltaW
        DEBUG_OFF('dyn-map {}'.format(i_track))
        DEBUG_OFF('dyn-map {}'.format(f_track))
        return f_track

    def adjust_energy(self, tkin):
        adjusted = DYNG(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted
class TestDynacGapMapping(unittest.TestCase):
    def test_DYNG_REF_mapping(self):
        print('----------------------------------test_DYNG_REF_mapping')
        label    = 'dyn_gap_test'
        phisoll  = MATH.radians(-25.)
        gap      = 0.044
        freq     = 800e6
        tkin     = 50.
        particle = Proton(tkin)
        dWf      = 1
        fname    = 'Superfish/SF_WDK2g44.TBL'
        gap_cm   = gap*100     # Watch out!
        EzPeak   = 10.0
        Ezdata   = SFdata(fname,EzPeak=EzPeak,gap=gap_cm)
        EzAvg    = Ezdata.EzAvg
        dyngap   = DYNG(label,EzAvg,phisoll,gap,freq,SFdata=Ezdata,particle=particle,dWf=dWf)
        i_track  = NP.array([0.,0.,0.,0.,0.,0.,tkin,1,0,1])
        f_track  = dyngap.dynac_map(i_track)
if __name__ == '__main__':
    unittest.main()
