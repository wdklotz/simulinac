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
import IGap
import unittest
import numpy as NP
from copy import copy
from math import sqrt,degrees,cos,sin,pi
from setutil import I0,I1,WConverter,Proton,wrapRED,PARAMS,FLAGS,Twiss,OutOfRadialBoundEx
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil import DEBUG_ON,DEBUG_OFF,log_what_in_interval
from separatrix import w2phi

twopi = 2.*pi

class Base_G(IGap.IGap):
    """ Base RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self):
        pass

    def configure(self,**kwargs):
        # self.length       = 0. # 0. because it's a kick
        self.dWf          = FLAGS['dWf']
        self.mapping      = 'base'        # map model
        self.kwargs       = kwargs
        self.label        = 'BM' 

        self.EzPeak    = kwargs.get('EzPeak',None)
        self.phisoll   = kwargs.get('phisoll',None)
        self.cavlen    = kwargs.get('cavlen',None)
        self.freq      = kwargs.get('freq',None)
        self.particle  = kwargs.get('particle',Proton(50.))
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)

        self.lamb      = PARAMS['clight']/self.freq
        self.gap       = self.cavlen*0.57   # 57%: a best guess ?
        self.matrix    = None
        self.ttf       = None
        self.deltaW    = None
        self.particlef = None
        self.master    = None

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return self.base_map(i_track)
        # return self.base_map_2(i_track)

    def toString(self):
        return f'{self.mapping} mapping in: Base_M.base_map()'

    def isAccelerating(self):
        return True

    def waccept(self):
        """ 
        Calculate longitudinal acceptance, i.e. phase space ellipse parameters: T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0*Dphi0 = ellipse_area/pi
        """
        # key-word parameters
        twiss_w_i = PARAMS['twiss_w_i']    # beta, alpha, gamma, eittance @ entrance (callable oject!)
        Dphi0_i   = PARAMS['Dphi0_i']      # Delta-phi @ entrance [rad]
        DT2T_i    = PARAMS['DT2T_i']       # Delta-T/T @ entrance

        # instance members
        Ez0       = self.EzPeak
        ttf       = self.ttf
        phisoll   = self.phisoll         # [rad]
        lamb      = self.lamb            # [m]
        freq      = self.freq            # [Hz]
        particle  = self.particle

        # calculated variables
        E0T       = Ez0*ttf              # [MV/m]
        m0c2      = particle.e0          # [MeV]
        gb        = particle.gamma_beta
        beta      = particle.beta
        gamma     = particle.gamma
        tkin      = particle.tkin

        # converter for this object 
        conv = WConverter(tkin,freq)

        try:
            # LARGE amplitude oscillations (T.Wangler pp. 175 6.28). w = Dgamma = DW/m0c2 normalized energy spread """
            # DEBUG_OFF(f'w2phi {(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll)}')                                                                                                                                                              
            w0large = sqrt(w2phi(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll))
        except ValueError as ex:
            exception = ex
            w0large = -1
        try:
            # SMALL amplitude oscillations separatrix (T.Wangler pp.185) """
            w0small = sqrt(2.*E0T*gb**3*lamb*phisoll**2*sin(-phisoll)/(pi*m0c2))
        except ValueError as ex:
            w0small = -1

        if w0large != -1: 
            wmax = w0large
        elif w0large == -1 and w0small != -1:
            wmax = w0small
        else:
            raise(UserWarning(wrapRED(f'{ex} reason: ttf={rf_gap.ttf}, E0T={E0T}')))
            sys.exit(1)

        # Dp/p max on separatrix
        Dp2pmax = conv.wToDp2p(wmax) 

        try:
            #  convert T.Wangler units {Dphi,w} to {z,dp/p} units with entrance parameters
            (beta_wsp,dummy,dummy,emit_wsp) = twiss_w_i()
            w0_i = (gamma-1.)*DT2T_i
            z0_i,Dp2p0_i,emitz_i,betaz_i = conv.wtoz((Dphi0_i,w0_i,emit_wsp,beta_wsp))
            alfaz_i = 0.
        except RuntimeError as ex:
            print(wrapRED(ex))
            sys.exit()

        # omega sync for this node
        omgl_0 = sqrt(E0T*lamb*sin(-phisoll)/(twopi*m0c2*gamma**3*beta))*twopi*freq   # [Hz]

        # phase acceptance (REMARK: phase limits are not dependent on Dp/p aka w)
        phi_2=2.*phisoll
        phi_1=-phisoll

        res =  dict (
                emitw_i         = emit_wsp,     # emittance {Dphi,w} units [rad,1]
                z0_i            = z0_i,         # ellipse z-axe crossing (1/2 axis) [m]
                Dp2p0_i         = Dp2p0_i,      # ellipse dp/p-axe crossing (1/2 axis)
                twiss_z_i       = Twiss(betaz_i, alfaz_i, emitz_i), # cavity twiss parameters
                DWmax           = wmax*m0c2,    # max delta-W on separatrix [MeV]
                Dp2pmax         = Dp2pmax,      # Dp/p max on separatrix [1]
                phaseacc        = (conv,phi_2,phisoll,phi_1), # phase acceptance [rad]
                omgl_0          = omgl_0,       # synchrotron oscillation [Hz]
                wmax            = wmax,         # w max on separatrix [1] (large amp. oscillations)
                zmax            = conv.DphiToz(-phisoll) # z max on separatrix [m] (large amp. oscillations -- Wrangler's approximation (pp.178) is good up to -58deg)
                )
        return res

    def register_mapper(self,master):
        master.register_mapping(self)
        pass

    def accept_register(self,master):
        self.master = master
        pass

    def adjust_energy(self, tkin):
        self.particle  = Proton(tkin)
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.deltaW    = self.EzPeak * self.ttf * self.gap * cos(self.phisoll)
        self.particlef = Proton(tkin + self.deltaW)
        self.T3D_matrix()
        pass

    def base_map(self, i_track):
        """ Neue überarbeitete map Version vom 17.02.2025 (wdk)
            Mapping in Base RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy Sollteilchen
        S        = i_track[SKOO]       # [8] position gap

        particleIs = self.particle   # Sollteilchen energy In
        m0c2        = particleIs.e0
        betaIs      = particleIs.beta
        gbIs        = particleIs.gamma_beta
        tkinIs      = particleIs.tkin
        freq        = self.freq
        lamb        = self.lamb
        phisoll     = self.phisoll
        L           = self.gap
        ttf         = self.ttf
        qE0LT       = self.EzPeak * L * ttf
        DtkinS      = self.deltaW
        tkinSo      = tkinIs + DtkinS  # Sollteilchen energy Out
        
        max_r  = 0.05              # max radial excursion [m]
        r      = sqrt(x**2+y**2)   # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbIs)
        i0     = I0(Kr)            # bessel function I0
        i1     = I1(Kr)            # bessel function I1

        # Teilchen
        converter   = WConverter(tkinIs,freq)
        phiI        = converter.zToDphi(z) + phisoll    # Phase Teilchen In
        Dtkin       = qE0LT*i0*cos(phiI)                # Teilchen energy kick (Shislo 4.2.3)
        tkinI       = converter.Dp2pToDW(zp) + tkinIs   # Teilchen energy In 
        tkinO       = tkinI + Dtkin                     # Teilchen energy Out   
        DDtkin      = tkinO - tkinSo                    # Differenz von Teichen und Sollteilchen Out

        particleOs  = Proton(tkinSo)
        betaOs      = particleOs.beta
        gbOs        = particleOs.gamma_beta

        zO          = betaOs/betaIs*z                     # z A.Shishlo/J.Holmes (4.2.5) 
        zpO         = converter.DWToDp2p(DDtkin)         # DDtkin --> dp/p

        # 17.02.2025 wdk: verbessertes Abfangen wenn lim(r)->0
        i12r = i1/r if r > 1.e-6 else 0.5*twopi/lamb/gbIs
        factor = qE0LT/(m0c2*gbIs*gbOs)*i12r*sin(phiI)  # common factor
        gbi2gbo = gbIs/gbOs
        xp  = gbi2gbo*xp - factor*x   # Formel 4.2.6 A.Shishlo/J.Holmes
        yp  = gbi2gbo*yp - factor*x

        f_track = NP.array([x, xp, y, yp, zO, zpO, T+DtkinS, 1., S, 1.])

        # S1 = 49.750    # from
        # S2 = 49.770    # to 
        # log_what_in_interval(S,(S1,S2),f'Base_M.base_map: f_track: {f_track}\n')

        return f_track

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

def ttf(lamb, gap, beta):
    """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
    x = gap/(beta*lamb)
    res =NP.sinc(x)
    return res

class TestBaseMapping(unittest.TestCase):
    def test_BASE_M(self):
        """ testing the base mapping for acceleration """
        print(wrapRED('------------------ test_BASE_M'))
        gap_parameter = dict(
            EzPeak    = 1,
            phisoll   = radians(-30.),
            cavlen    = 0.44,
            freq      = 750e6,
        )
        bmap = BASE_M()    # create object instance
        bmap.configure(**gap_parameter)
        print(bmap.toString())

if __name__ == '__main__':
    unittest.main()

