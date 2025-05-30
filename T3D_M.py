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
from setutil    import PARAMS,FLAGS,Proton,MDIM,DEBUG_ON,DEBUG_OFF,mxprnt,wrapRED
from setutil    import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil    import WConverter,Twiss,I0
from math       import pi,radians,sin,cos,sqrt
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


class T3D_G(IGap.IGap):
    """ Trace 3D RF Gap-Model """
    def __init__(self):
        self.label        = 'T3D_G'
        self.master       = None

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
    def ttf(self):           return self.master.ttf             # ttf
    @ttf.setter
    def ttf(self,v):                self.master.ttf = v

    def accelerating(self):     return self.master.accelerating
    def adjust_energy(self,tkin):
        self.T3D_matrix()
        pass
    def configure(self,**kwargs):
        # copies from master
        self.aperture  = kwargs.get('aperture')
        self.cavlen    = kwargs.get('cavlen')
        self.EzPeak    = kwargs.get('EzPeak')
        self.freq      = kwargs.get('freq')
        self.gap       = kwargs.get('gap')
        self.phisoll   = kwargs.get('phisoll')
        self.sec       = kwargs.get('sec')

        self.lamb      = kwargs['lamb']
        self.omega     = kwargs['omega'] 
        self.E0L       = None
        self.qE0LT     = None
    def map(self,i_track):      return NP.dot(self.matrix,i_track)
    def register(self,master):
        self.master = master
        pass
    def toString(self):         return mxprnt(self.matrix,'4g')
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
            #DEBUG_OFF(f'w2phi {(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll)}')                                                                                                                                                              
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
        m              = NP.eye(MDIM,MDIM)
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta,self.aperture)
        self.deltaW    = self.E0L*self.ttf*cos(self.phisoll)
        self.particlef = Proton(self.particle.tkin+self.deltaW)
        Wavg    = self.particle.tkin+self.deltaW/2.   # average tkin
        pavg    = Proton(Wavg)
        bavg    = pavg.beta
        gavg    = pavg.gamma
        m0c2    = pavg.e0
        kz      = twopi*self.E0L*self.ttf*sin(self.phisoll)/(m0c2*bavg*bavg*self.lamb)
        ky      = kx = -0.5*kz/(gavg*gavg)
        bgi     = self.particle.gamma_beta
        bgf     = self.particlef.gamma_beta
        bgi2bgf = bgi/bgf
        m       = NP.eye(MDIM,MDIM)
        m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
        m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
        m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
        m[EKOO, DEKOO] = self.deltaW
        m[SKOO, DSKOO]  = 0.
        self.matrix = m
        pass
        return

class TestElementMethods(unittest.TestCase):
    def testT3D(self):
        """ testing the Trace 3D mapping for acceleration """
        print(wrapRED('--------------------------------------- T3D_G.T3D_matrix()'))
        gap_parameter = dict(
            aperture   = 0.01,
            cavlen     = 0.044,
            EzPeak    = 1.1,
            freq       = 750.e+6,
            gap        = 0.022,
            phisoll    = radians(-30.),
            sec        = 'test',
            SFdata     = None,
        )
        FLAGS['mapping'] = 't3d'
        instance = ELM.RFG('RFG',PARAMS['injection_energy'])
        instance.register(T3D_G())
        instance.configure(**gap_parameter)
        instance.adjust_energy(PARAMS['injection_energy'])

        tkin = instance.particle.tkin
        print(f'transfer matrix T3D for {tkin} MeV (default)')
        print(instance.mapper.toString())

        print(wrapRED('------------------------------------------ T3D_G.waccept()'))
        DT2T       = 2e-3  # DT/T
        Dphi0      = 5.    # Dphi/phi [deg]
        E0         = PARAMS['proton_mass']
        w0         = tkin/E0*DT2T # Wrangler's definition of w (pp.176)
        emit       = Dphi0*w0     # emittance  in {Dphi,w}-space
        betaw      = emit/w0**2   # twiss-beta in {Dphi,w}-space
        PARAMS['twiss_w_i'] = Twiss(beta=betaw,alfa=0.,epsi=emit)
        PARAMS['Dphi0_i']   = Dphi0
        PARAMS['DT2T_i']    = DT2T
        pprint.pp(instance.waccept())
        print('twiss_z_i')
        pprint.pp(vars(instance.waccept()['twiss_z_i']))

if __name__ == '__main__':
    import elements as ELM
    import pprint
    unittest.main()
