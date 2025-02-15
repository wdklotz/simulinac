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
import IGap
import unittest
import numpy as NP
from setutil import PARAMS,FLAGS,Proton,MDIM,DEBUG_ON,DEBUG_OFF,mxprnt,wrapRED
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil import WConverter,Twiss,dictprnt
from math import pi,radians,sin,cos,sqrt
from separatrix import w2phi

twopi = 2.*pi

def ttf(lamb, gap, beta):
    """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
    x = gap/(beta*lamb)
    res =NP.sinc(x)
    return res

class T3D_G(IGap.IGap):
    """ Trace 3D RF Gap-Model """
    def __init__(self):
        pass

    def configure(self,**kwargs):
        self.length       = 0. # 0. because it's a kick
        self.dWf          = FLAGS['dWf']                 # dWf=1 with acceleration =0 else

        self.kwargs    = kwargs
        self.label     = kwargs.get('label','3D')  
        self.EzPeak    = kwargs.get('EzPeak',None)*self.dWf # [MV/m]
        self.phisoll   = kwargs.get('phisoll',None)         # [radians] soll phase
        self.gap       = kwargs.get('gap',None)             # [m] rf-gap
        self.cavlen    = kwargs.get('cavlen',None)          # [m] cavity length
        self.freq      = kwargs.get('freq',None)            # [Hz]  RF frequenz
        self.particle  = kwargs.get('particle',Proton(50.))
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)
        self.mapping   = kwargs.get('mapping','t23d')        # map model

        self.omega     = twopi*self.freq
        self.lamb      = PARAMS['clight']/self.freq
        self.ttf       = None
        self.E0L       = None
        self.qE0LT     = None
        self.deltaW    = None
        self.particlef = None
        self.matrix    = None
        self.T3D_matrix()
        pass

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return NP.dot(self.matrix,i_track)
    
    def toString(self):
        return mxprnt(self.matrix,'4g')

    def isAccelerating(self):
        return True

    def adjust_energy(self,tkin):
        self.particle = Proton(tkin)
        self.T3D_matrix()
        pass

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

    def T3D_matrix(self):
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m              = NP.eye(MDIM,MDIM)
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.qE0LT     = self.E0L*self.ttf
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
        return

class TestElementMethods(unittest.TestCase):
    def testT3D1(self):
        """ testing the Trace 3D mapping for acceleration """
        print(wrapRED('------------------ testT3D1'))
        gap_parameter = dict(
            EzPeak    = 1,
            phisoll   = radians(-30.),
            gap       = 0.022,
            cavlen    = 0.44,
            freq      = 750e6,
        )
        t3dg = T3D_G()    # create object instance
        t3dg.configure(**gap_parameter)
        print(f'transfer matrix T3D for {t3dg.particle.tkin} MeV (default)')
        print(t3dg.toString())

        tkin = 100.
        print(f'transfer matrix T3D for {tkin} MeV')
        t3dg.adjust_energy(tkin)
        print(t3dg.toString())
        #=============================================================
        tkin       = 10
        print(f'\n test waccept(...) for {tkin} MeV protons')
        DT2T       = 2e-3  # DT/T
        Dphi0      = 5.    # Dphi/phi [deg]
        E0         = PARAMS['proton_mass']
        w0         = tkin/E0*DT2T # Wrangler's definition of w (pp.176)
        emit       = Dphi0*w0     # emittance  in {Dphi,w}-space
        betaw      = emit/w0**2   # twiss-beta in {Dphi,w}-space
        PARAMS['twiss_w_i'] = Twiss(beta=betaw,alfa=0.,epsi=emit)
        PARAMS['Dphi0_i']   = Dphi0
        PARAMS['DT2T_i']    = DT2T
        res = t3dg.waccept()
        dictprnt(what=res,text='waccept',njust=25)

if __name__ == '__main__':
    unittest.main()
