#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='v11.0.3'
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
from math       import sqrt,radians,cos,sin,pi
from setutil    import I0,I1,WConverter,Proton,wrapRED,PARAMS,FLAGS,Twiss,OutOfRadialBoundEx
from setutil    import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil    import DEBUG_ON,DEBUG_OFF,mxprnt,log_what_in_interval
from separatrix import w2phi

twopi = 2.*pi
def ttf(lamb, gap, beta, aperture):
    """ WRANGLER: Transit-Time-Factor Models, pp. 44 (2.43) """
    x   = gap/(beta*lamb)
    res = NP.sinc(x)
    gamma = 1/sqrt(1-beta**2)
    Ka = twopi/(lamb*gamma*beta)*aperture
    res = res/I0(Ka)
    return res

class BASE_G(IGap.IGap):
    """ BASE RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self):
        self.master       = None
        self.label        = 'BASE_G'

    # mutable properties shared with master
    @property
    def deltaW(self):        return self.master.deltaW          # deltaW
    @deltaW.setter
    def deltaW(self,v):             self.master.deltaW = v
    @property
    def matrix(self):        return self.master.matrix          # matrix
    @matrix.setter
    def matrix(self,v):             self.master.matrix = v
    @property
    def particle(self):      return self.master.particle        # particle
    @property
    def particlef(self):     return self.master.particlef       # particlef
    @particlef.setter
    def particlef(self,v):          self.master.particlef = v
    @property
    def ttf(self):           return self.master.ttf              # ttf
    @ttf.setter
    def ttf(self,v):                self.master.ttf = v

    def accelerating(self):    return self.master.accelerating
    def adjust_energy(self, tkin):
        self.ttf = ttf(self.lamb, self.gap, self.particle.beta, self.aperture)
        self.T3D_matrix()
        pass
    def configure(self,**kwargs):
        self.aperture  = kwargs.get('aperture')
        self.cavlen    = kwargs.get('cavlen')
        self.EzPeak    = kwargs.get('EzPeak')
        self.freq      = kwargs.get('freq')
        self.gap       = kwargs.get('gap')
        self.HE_Gap    = kwargs.get('HE_Gap')
        self.phisoll   = kwargs.get('phisoll')
        self.sec       = kwargs.get('sec')
        self.SFdata    = kwargs.get('SFdata')

        self.lamb      = kwargs['lamb']
        self.omega     = kwargs['omega'] 
    def map(self,i_track):     return self.base_map(i_track)
    def register(self,master):
        self.master = master
        pass
    def toString(self):        return mxprnt(self.matrix,'4g')
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
            exception = ex
            w0small = -1

        if w0large != -1: 
            wmax = w0large
        elif w0large == -1 and w0small != -1:
            wmax = w0small
        else:
            raise(UserWarning(wrapRED(f'{exception} reason: ttf={self.ttf}, E0T={E0T}')))
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
    def T3D_matrix(self):
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m         = NP.eye(MDIM,MDIM)
        E0L       = self.EzPeak*self.gap
        deltaW    = self.deltaW = E0L*self.ttf*cos(self.phisoll)
        Wavg      = self.particle.tkin+self.deltaW/2.   # average tkin
        pavg      = Proton(Wavg)
        bavg      = pavg.beta
        gavg      = pavg.gamma
        m0c2      = pavg.e0
        kz        = twopi*E0L*self.ttf*sin(self.phisoll)/(m0c2*bavg*bavg*self.lamb)
        ky        = kx = -0.5*kz/(gavg*gavg)
        bgi       = self.particle.gamma_beta
        particlef = Proton(self.particle.tkin + deltaW)
        bgf       = particlef.gamma_beta
        bgi2bgf   = bgi/bgf
        m         = NP.eye(MDIM,MDIM)
        m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
        m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
        m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf
        m[EKOO, DEKOO] = deltaW
        m[SKOO, DSKOO]  = 0.
        self.matrix = m
        return

    def base_map(self, i_track):
        """ Neue überarbeitete map Version vom 17.02.2025 (wdk)
            Mapping in BASE RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kin energy Soll
        S        = i_track[SKOO]       # [8] position gap

        particleIs  = self.particle   # Soll energy In
        m0c2        = particleIs.e0
        betaIs      = particleIs.beta
        gbIs        = particleIs.gamma_beta
        tkinIs      = particleIs.tkin
        freq        = self.freq
        lamb        = self.lamb
        phisoll     = self.phisoll
        qE0LT       = self.EzPeak * self.gap * self.ttf
        DtkinS      = self.deltaW
        tkinSo      = tkinIs + DtkinS  # Soll energy Out
        
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
        DDtkin      = tkinO - tkinSo                    # Differenz von Teichen und Soll Out

        particleOs  = Proton(tkinSo)
        betaOs      = particleOs.beta
        gbOs        = particleOs.gamma_beta

        zO          = betaOs/betaIs*z                     # z (Shishlo 4.2.5) 
        zpO         = converter.DWToDp2p(DDtkin)          # DDtkin --> dp/p

        # 17.02.2025 wdk: verbessertes Abfangen wenn lim(r)->0
        i12r = i1/r if r > 1.e-6 else pi/(lamb*gbIs)
        factor = qE0LT/(m0c2*gbIs*gbOs)*i12r*sin(phiI)  # common factor
        gbi2gbo = gbIs/gbOs
        xp  = gbi2gbo*xp - factor*x   # Formel 4.2.6 A.Shishlo/J.Holmes
        yp  = gbi2gbo*yp - factor*y

        f_track = NP.array([x, xp, y, yp, zO, zpO, T+DtkinS, 1., S, 1.])

        S1 = 48    # from
        S2 = 51    # to 
        debug = DEBUG_OFF
        if debug == DEBUG_ON:
            log_what_in_interval(S,(S1,S2),f'BASE_M.base_map: f_track: {f_track}\n')

        return f_track

class TestBASEMapping(unittest.TestCase):
    def test_BASE_G(self):
        """ testing the base mapping for acceleration """
        print(wrapRED('---------------------------------------------------------- test_BASE_G'))
        gap_parameter = dict(
            EzPeak    = 2,
            phisoll   = radians(-30.),
            gap       = 0.044,
            aperture  = 0.010,
            freq      = 750e6,
        )
        FLAGS['mapping'] = 'base'
        instance = ELM.RFG('RFG')
        instance.register(BASE_G())
        instance.configure(**gap_parameter)
        instance.adjust_energy(PARAMS['injection_energy'])

        print(instance.mapper.toString())
        print(instance.mapper.__dict__)
        print()
        x  = 1e-3
        xp = 1e-3
        y  = 1e-3
        yp = 1e-3
        z  = 10e-3
        zp = 1e-3
        T  = 50
        S  = 99
        intrack=NP.array([x, xp, y, yp, z, zp, T, 1, S, 1])
        outrack1=instance.map(intrack)
        outrack2=NP.dot(instance.matrix,intrack)
        prec=9
        print(f'map: {NP.array_repr(outrack1,precision=prec,suppress_small=True,max_line_width=180)}')
        print(f't3d: {NP.array_repr(outrack2,precision=prec,suppress_small=True,max_line_width=180)}')
        print(f' -   {NP.array_repr(outrack1-outrack2,precision=prec,suppress_small=True,max_line_width=180)}')

if __name__ == '__main__':
    import elements as ELM
    unittest.main()

