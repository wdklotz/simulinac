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
from functools import partial
import warnings
from collections import namedtuple

import elements as ELM
from setutil import DEB, arrprnt, PARAMS, tblprnt, Ktp, WConverter
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Ipoly, Polyval

DEBUG_OFF = DEB.get('OFF')
DEBUG_ON  = DEB.get('ON')

twopi = 2.*MATH.pi

""" Numerical computations in an accelerating gap or in a cavity; E.TANKE and S.VALERO 3-Oct-2016 """

class StepFactory(object):
    """ StepFactory """
    def __init__(self,gap,SFdata):
        self.zl = -gap*100./2. # [cm]
        self.zr = -self.zl
        polyvals = []
        for poly in SFdata.EzPoly:
            zil = poly.zl
            zir = poly.zr
            if zil < self.zl or zir > self.zr: 
                continue
            else:
                polyvals.append((poly.zl*1.e-2,poly.zr*1.e-2))    # [m]
        polyvals = tuple(polyvals) 
        self.h = polyvals[0][1] - polyvals[0][0]
        z_parts = []
        for interval in polyvals:
            z0 = interval[0]
            z1 = z0 + self.h/4.
            z2 = z0 + self.h/2.
            z3 = z0 + (3*self.h)/4.
            z4 = z0 + self.h
            z_parts.append((z0,z1,z2,z3,z4))
        self.z_steps = NP.array(z_parts)
        return
    def zArray(self,n):
        return self.z_steps[n]
    def tArray(self,t0,betac):
        t1 = t0 + self.h/(4*betac)
        t2 = t0 + self.h/(2*betac)
        t3 = t0 +(3*self.h)/(4*betac)
        t4 = t0 + self.h/betac
        return NP.array([t0,t1,t2,t3,t4])
    def nsteps(self):
        return len(self.z_steps)
    def nparts(self):
        return 5
    def steplen(self):
        return self.h
class DYNG(ELM.RFG):
    """ DYNAC's RF-gap model """
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
        super().__init__(label, EzAvg, phisoll, gap, freq, particle, position, aperture, dWf, mapping='dyn')
        if SFdata == None:
            raise RuntimeError('DYNG: missing E(z) table - STOP')
            sys.exit(1)
        else:
            self.SFdata    = SFdata
            self.map       = self.dynac_map   # OXAL's specific mapping method
            self.polies    = self.poly_slices(self.gap,self.SFdata)
            self.particlef = Proton(particle.tkin + self.deltaW)
            # TODO self.deltaW    = self.matrix[ Ktp.T,Ktp.dT]
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
            sqrt_g2m1      = sqrt(g2msqrt)
            sqrt_sqrt_g2m1 = sqrt(sqrt_g2m1)
            if not back:
                # r=(x,y) rp=(xp,yp)
                X = sqrt_sqrt_g2m1*x
                Y = sqrt_sqrt_g2m1*y
                Xp = (sqrt_sqrt_g2m1*xp+0.5*X*gamma/g2m1)
                Yp = (sqrt_sqrt_g2m1*yp+0.5*y*gamma/g2m1)
                return (X,XP,Y,YP)
            else:
                # r=(X,Y) rp=(Xp,Yp)
                X = x; XP=xp; Y=y; YP=yp
                x = X/sqrt_sqrt_g2m1
                y = Y/sqrt_sqrt_g2m1
                xp = (XP-0.5*gamma*X/g2m1)/sqrt_sqrt_g2m1)
                yp = (YP-0.5*gamma*Y/g2m1)/sqrt_sqrt_g2m1)
                return (x,xp,y,yp)
    #TODO: use non_const beta, i.e. gamma[i], i=1,4,1 in inegrals
    #TODO: find position for multiplication with dWf-flag
    #TODO: use equivalent E-field instead of SFdata
    def dynac_map(self,i_track):
        """ Mapping from (i) to (f) """  
        # stpfac = self.stpfac
        c      = PARAMS['clight']
        h      = stpfac.steplen()  
        m0c2   = PARAMS['proton_mass']
        m0c3   = m0c2*c
        freq   = self.freq
        omega  = self.omega
        # initial values for loop over polies
        p           = copy(self.particle)
        Ws_in       = p.tkin
        phis_in     = self.phisoll
        f_track     = copy(i_track)
        self.deltaW = -Ws_in
        for poly in self.polies:    
            """ Map through a single poly interval 
            For each poly interval a numerical step-by-step method based on the 5 point Boole's will be applied. The step size in the azimuthally direction z is divided in 4 parts equivalent 
            lengths."""     
            # Azimutally positions z0,z1,z2,z3,z4; step size h
            h  = poly.dz*1.e-2     # stepsize [m]
            z0 = poly.zl*1.e-2     # [m]
            z1 = z0 + h/4.
            z2 = z0 + h/2.
            z3 = z0 = 3.*h/4.
            z4 = z0 + h
            
            gbs_in = p.gamma_beta
            gbs03 = gbs_in**3
            K1 = omega**2/4/c**2/gbs03
            x         = f_track[XKOO]       # [0]
            xp        = f_track[XPKOO]      # [1]
            y         = f_track[YKOO]       # [2]
            yp        = f_track[YPKOO]      # [3]
            gamma_in   = p.gamma`
            Rx_in,Rxp_in,Ry_in,Ryp_in = self.Picht(gamma0,x,xp,y,yp,back=False)
            betas_in = p.beta
            bc       = betas_in*c
            t0 = z0/bc
            t1 = t0 + h/(4*bc)
            t2 = t0 + h/(2*bc)
            t3 = t0 + 3.*h/(4*bc)
            t4 = t0 + h/bc

            Ez0 = self.SFdata.Ez0t(z0*100.,t0,omega,phis_in)
            Ez1 = self.SFdata.Ez0t(z1*100.,t1,omega,phis_in)
            Ez2 = self.SFdata.Ez0t(z2*100.,t2,omega,phis_in)
            Ez3 = self.SFdata.Ez0t(z3*100.,t3,omega,phis_in)
            Ez4 = self.SFdata.Ez0t(z4*100.,t4,omega,phis_in)
            I1  = h/90.*(7.*Ez0 + 32.*Ez1 + 12.*Ez2 + 32.*Ez3 + 7.*Ez4)
            I2  = h**2/90.*(8.*Ez1 + 6.*Ez2 + 24.* Ez3 + 7.*Ez4)
            R02   = Rx_in**2 + Ry_in**2
            R0R0p = Rx_in*Rxp_in + Ry_in*Ryp_in
            """ Energy """
            Dgamma = 1/m0c2*((1+R02*K1)*I1 + R0R0p*K1*I2)
            gamma_out = gamma0 + Dgamma

            I3 = h**2/90./gbs03*(8.*Ez1 + 6.*Ez2 + 24.*Ez3 + 7.*Ez4)
            I4 = h**3/90./gbs03*(2.*Ez1 + 3.*Ez2 + 18.*Ez3 + 7.*Ez4)
            """ Phase """
            Dt = 1/m0c3*((1+R02*K1)*I3 + R0R0p*K1*I4)
            t4_out = t4 + Dt
            """ transverse x,xp,y,yp """
            g2m1 = gamma_in-1
            G1 = 1./2/m0c3*sqrt(1./g2m1**3)
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

            J3 = H**3/90./G1*((2.*Ez1p + 3.*Ezp2 + 18.*Ezp3 + 7.*Ezp4))

            DRx_out = Rx_in*J2 + Rxp_in*J3
            Rx_out = Rx_in + DRx_out + h*Rxp_in
            DRy_out = Ry_in*J2 + Ryp_in*J3
            Ry_out = Ry_in + DRy_out + h*Ryp_in

            x,xp,y,yp = self.Picht(gamma_out,Rx_out,Rxp_out,Ry_out,Ryp_out,back=True)

            """ Reset der loop Variablen  ??????????????????????? """


            # TODO ------------------------------------------------------
        f_track = NP.array([ x, xp, y, yp, z, zp, T, 1., S, 1.])
        DEBUG_OFF('dyn-map {}'.format(i_track))
        DEBUG_OFF('dyn-map {}'.format(f_track))
        return f_track
def test0():
    DEBUG_TEST0 = DEBUG_ON
    import elements as ELM
    print('-----------------------------------TEST 0----------------')
    print('test _DYN_Gslice:slice_map()...')
    input_file='SF_WDK2g44.TBL'
    Ezpeak = 1.4
    SF_tab = SFdata(input_file,Ezpeak)
    x  = 1.0e-2
    xp = 1.0e-3
    y  = 1.0e-2
    yp = 1.0e-3
    z  = 1.0e-3
    zp = 1.0e-3
    
    i_track1 = NP.array([ 0,  0, 0,  0, 0,  0, PARAMS['sollteilchen'].tkin, 1., 0., 1.])
    i_track2 = NP.array([ x, xp, y, yp, z, zp, PARAMS['sollteilchen'].tkin, 1., 0., 1.])
    dyng = ELM.RFC(gap=0.048, SFdata=SF_tab, mapping='dyn')
    
    f_track = dyng.soll_map(i_track1)
    DEBUG_TEST0('dyn-soll:i_track:\n{}'.format(str(i_track1)))
    DEBUG_TEST0('dyn-soll:f_track:\n{}'.format(str(f_track)))

    DEBUG_TEST0('dyn-map:i_track:\n{}'.format(str(i_track2)))
    tracks = []
    for i in range(10):
        f_track = dyng.map(i_track2)
        tracks.append(f_track)
        i_track2 = f_track
    DEBUG_TEST0('dyn-map:f_track:')
    for track in tracks:
        DEBUG_TEST0('{:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} {:+10.5f} '.format(track[0],track[1],track[2],track[3],track[4],track[5],track[6],track[7],track[8],track[9]))
if __name__ == '__main__':
    test0()
