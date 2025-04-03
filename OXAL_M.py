#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='11.0.2.4'
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
from math import sin,cos,tan,sqrt,pi,degrees,radians
from setutil import FLAGS,PARAMS,Ktp,MDIM,Proton,DEBUG_ON,DEBUG_OFF
from setutil import wrapRED,mxprnt,Twiss,WConverter,dictprnt
from Ez0 import SFdata
from separatrix import w2phi
import Ez0 as EZ
# from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM


twopi = 2*pi
pihalf = pi/2
counter_of_polies = 0
trigger_poly_number = 12
twothird = 2./3.

def cot(x):
    return -tan(x+pihalf)

class OXAL_G(IGap.IGap):
    """ OpenXAL RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self):
        self.master      = None
        self.label       = 'OXAL_G'

    def configure(self,**kwargs):
        self.aperture  = kwargs.get('aperture')
        self.cavlen    = kwargs.get('cavlen')
        self.EzPeak    = kwargs.get('EzPeak')
        self.freq      = kwargs.get('freq')
        self.phisoll   = kwargs.get('phisoll')
        self.sec       = kwargs.get('sec')
        self.SFdata    = kwargs.get('SFdata')

        self.lamb      = kwargs['lamb']
        self.omega     = kwargs['omega'] 
        self.polies    = self.poly_slices()

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

    def map(self,i_track):
        return NP.dot(self.matrix,i_track)
    def toString(self):
        return mxprnt(self.matrix,'4g')
    def isAccelerating(self):
        return True
    def adjust_energy(self, tkin):
        self.OXAL_matrix()
        pass
    def waccept(self,**kwargs):
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
            # DEBUG_ON(f'w2phi {(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll)}')                                                                                                                                                              
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
    def register(self,master):
        self.master = master
        pass
    def OXAL_matrix(self):
        polies   = self.polies
        c        = PARAMS['clight']
        m0c2     = self.particle.m0c2
        m0c3     = m0c2*c
        omega    = self.omega
        ttf      = 0.
        matrix   = NP.eye(MDIM,MDIM)

        # initialise loop variables
        Ts       = self.particle.tkin                 # T
        phis     = self.phisoll                       # phi

        for poly in polies:   # each poly is a slice of the full Ez-dist
            # IN variables
            gammas_in  = 1. + Ts/m0c2              # gamma 
            gbs_in     = sqrt(gammas_in**2-1.)     # (gamma*beta)
            betas_in   = gbs_in/gammas_in          # beta
            gbs3_in    = gbs_in**3                  # (gamma*beta)**3
            g3b2s_in   = gammas_in**3*betas_in**2   # gamma**3*beta**2
            g2s_in     = gammas_in**2               # gamma**2

            ks     = omega/(c*betas_in)    # omega/(beta*c)
            qV0    = self.V0(poly)         # [MV]
            Tks    = self.T(poly,ks)
            Sks    = self.S(poly,ks)
            Tpks   = self.Tp(poly,ks)
            Spks   = self.Sp(poly,ks)
            Tppks  = self.Tpp(poly,ks)  
            Sppks  = self.Spp(poly,ks) 

            sphis  = sin(phis)
            cphis  = cos(phis)

            # kin.energy increase ref-particle
            dTs       = qV0*(Tks*cphis - Sks*sphis)                   # Shishlo 4.6.1
            # phase increase ref-particle
            dPhis = qV0*omega/m0c3/gbs3_in*(Tpks*sphis + Spks*cphis)  # Shishlo 4.6.2
            Ts_out    = Ts + dTs
            phis_out  = phis + dPhis
            ttf = ttf + qV0

            # OUT variables
            gammas_out  = 1. + Ts_out/m0c2            # gamma 
            gbs_out     = sqrt(gammas_out**2-1.)      # (gamma*beta)
            betas_out   = gbs_out/gammas_out          # beta
            g3b2s_out   = gammas_out**3*betas_out**2  # gamma-s**3*beta-s**2
            g2s_out     = gammas_out**2               # gamma**2
            #======================================================================================="""        
            # OXAL-matrix 
            # (4.6.11) in Shishlo's paper:
            factor  = qV0*betas_out/m0c2/gbs3_in
            factor1 = qV0*omega/m0c3/betas_in/g3b2s_in
            r44 = betas_out/betas_in + factor*omega/c/betas_in*(Tpks*cphis - Spks*sphis)
            r45 = factor*(3*gammas_in**2*(Tpks*sphis + Spks*cphis) + omega/c/betas_in*(Tppks*sphis+Sppks*cphis))
            # (4.6.9) in Shishlo's paper:
            r54 = factor1*(Tks*sphis + Sks*cphis)
            r55 = (g3b2s_in/g3b2s_out - factor1*(Tpks*cphis - Spks*sphis))
            #======================================================================================="""        
            # {z, dP/P}: linear sub matrix
            # NOTE: Shishlo's Formeln sind fuer (z,dBeta/Beta)
            mx = NP.eye(MDIM,MDIM)
            mx[Ktp.z, Ktp.z] = r44;         mx[Ktp.z, Ktp.zp ] = r45/g2s_in  # g2s_in: apply conversion dBeta/Beta=gamma**(-2)*dP/P
            mx[Ktp.zp,Ktp.z] = r54*g2s_out; mx[Ktp.zp, Ktp.zp] = r55         # g2s_out: apply conversion dP/P=gamma**2*dBeta/Beta
            # {x,x'}: linear sub-matrix
            factor2 = qV0*omega/(2.*m0c3*gbs_out*gbs_in**2)
            mx[Ktp.xp,Ktp.x ] = -factor2 * (Tks*sphis + Sks*cphis)
            mx[Ktp.xp,Ktp.xp] = gbs_in/gbs_out
            # {y,y'}: linear sub-matrix
            mx[Ktp.yp,Ktp.y]  = mx[Ktp.xp, Ktp.x]
            mx[Ktp.yp,Ktp.yp] = mx[Ktp.xp, Ktp.xp]
            # energy and length increase
            mx[Ktp.T,Ktp.dT] = dTs
            mx[Ktp.S,Ktp.dS] = 0     # 0 length: oxal-gap is kick

            # left multiplication of slice-matrix with oxal-matrix
            matrix = NP.dot(mx,matrix)

            # refresh loop variables
            Ts    = Ts_out
            phis  = phis_out
        self.deltaW    = matrix[Ktp.T,Ktp.dT]
        self.ttf       = self.deltaW/ttf
        self.particlef = Proton(self.particle.tkin+self.deltaW)
        self.matrix    = matrix
        return
    def poly_slices(self):
        """Slice the RF cavity"""
        L = self.cavlen/2.
        slices = []
        zl = -L*100.   # [m] --> [cm]
        zr = -zl
        for poly in self.SFdata.polies:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slices.append(poly)
        DEBUG_OFF('slices',slices)
        return slices

    def V0(self, poly):      # A.Shishlo/J.Holmes (4.4.3)
        """ V0 A.Shishlo/J.Holmes (4.4.3) """
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = E0*(2*dz+twothird*b*dz**3)*1.e-2 # [MV]
        return v0 
    def T(self, poly, k):    # A.Shishlo/J.Holmes (4.4.6)
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+twothird*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz*cot(k*dz))
        t  = f1*f2
        DEBUG_OFF('TTF_G: (T,k) {}'.format((t,k)))
        return t
    def S(self, poly, k):    # A.Shishlo/J.Holmes (4.4.7)
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+twothird*b*dz**3))
        f2 = 1.-k*dz*cot(k*dz)
        s  = f1*f2
        DEBUG_OFF('TTF_G: (S,k) {}'.format((s,k)))
        return s
    def Tp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+twothird*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz*cot(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp
    def Sp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.9)
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        sp  = 2*a*sin(k*dz)/(k*(2*dz+twothird*b*dz**3))
        sp  = sp*(dz**2-2./k**2+dz*cot(k*dz)*2/k)
        sp  = sp*1.e-2     # [cm] --> [m]
        return sp
    def Tpp(self, poly, k):  # sympy calculation Tpp.ipynb
        """ 2nd derivative T''(k) """
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        I1 = 2*dz**2*sin(dz*k)/k + 4*dz*cos(dz*k)/k**2 - 4*sin(dz*k)/k**3
        I3 = -b*(-dz**4*sin(dz*k)/k - 4*dz**3*cos(dz*k)/k**2 + 12*dz**2*sin(dz*k)/k**3 + 24*dz*cos(dz*k)/k**4 - 24*sin(dz*k)/k**5) + b*(dz**4*sin(dz*k)/k + 4*dz**3*cos(dz*k)/k**2 - 12*dz**2*sin(dz*k)/k**3 - 24*dz*cos(dz*k)/k**4 + 24*sin(dz*k)/k**5)
        r = I1+I3
        return r
    def Spp(self, poly, k):  # sympy calculation Spp.ipynb
        """ 2nd derivative S''(k) """
        a   = poly.a
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        r = 2*dz**3*a*cos(dz*k)/k - 6*dz**2*a*sin(dz*k)/k**2 - 12*dz*a*cos(dz*k)/k**3 + 12*a*sin(dz*k)/k**4
        return r


class TestOxalEnergyMapping(unittest.TestCase):
    def test_OXAL(self):
        """ testing the OXAL mapping for acceleration """
        print(wrapRED('---------------------------------------- OXAL_G.OXAL_matrix()'))
        gap_parameter = dict(
            aperture   = 0.01,
            cavlen     = 0.044,
            EzPeak     = 1.1,
            freq       = 750.e+6,
            gap        = 0.022,
            phisoll    = radians(-30.),
            sec        = 'test',
            SFdata     = None,
        )
        FLAGS['mapping'] = 'oxal'
        fieldtab  = 'unittests/CAV-FLAT-R135-L32.TBL'
        L         = gap_parameter['cavlen']
        EzPeak    = gap_parameter['EzPeak']
        sfdata    = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=L*100.)   # scaled field distribution
        gap_parameter['SFdata'] = sfdata
        instance = ELM.RFG('RFG')
        instance.register(OXAL_G())
        instance.configure(**gap_parameter)
        instance.adjust_energy(PARAMS['injection_energy'])

        tkin = instance.particle.tkin
        print(f'transfer matrix OXAL for {tkin} MeV (default)')
        print(instance.mapper.toString())

        print(wrapRED('-------------------------------------------- OXAL_G.waccept()'))
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



""" altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug """
""" altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug """
""" altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug altes Zeug """
class KW():
    def OXAL_matrix(self):
        #====================================================================== sympy ==========
        #     """ k = omega/(c*beta), b = beta, g = gamma, om=omega
        #         Tk,Sk,Tpk,Spk,Tppk,Sppk = T(k),S(k),T'(k),S'(k),T''(k),S''(k) Fourierfaktoren
        #         Formeln mit sympy gerechnet Shishlo-OpenXal-1.ipynb
        #     """
        #     g2   = g**2      # gamma**2
        #     gb3  = (g*b)**3  # (gamma*beta)**3
        #     cphi = cos(phis)
        #     sphi = sin(phis)  

        #     # WoWi = qV0/g2*(     # Shishlo 4.6.1
        #     #   (Sk*cphi     - Tk*sphi)*g2*k*z        # ~z
        #     # - (Spk*sphi    - Tpk*cphi)*k*zp         # ~zp
        #     # - (Tk*cphi     - Sk*sphi)*g2
        #     # # - (Spk*cphi    - Tpk*sphi)*k**2*z*zp    # ~z*zp
        #     # )
        #     WoWi = qV0*(Sk*cphi*g2*k*z 
        #     - Sk*g2*sphi 
        #     # - Spk*cphi*k**2*z*zp 
        #     + Spk*k*sphi*zp 
        #     + Tk*cphi*g2 
        #     + Tk*g2*k*sphi*z 
        #     - Tpk*cphi*k*zp 
        #     # - Tpk*k**2*sphi*z*zp
        #     )/g2

        # PoPi = om*qV0/(b**3*g*g2**2*gb3*m0c3)*(b**3*g*g2*(   # Shishlo 4.6.2
        #   (Spk*cphi  + Tpk*sphi)*g2
        # + (Spk*sphi  - Tpk*cphi)*g2*k*z         # ~z
        # - (Sppk*cphi + Tppk*sphi)*k*zp          # ~zp
        # # - (Sppk*sphi - Tppk*cphi)*k**2*z*zp     # ~z*zp
        #   )
        # - 3*gb3*zp*( 
        #   (Spk*cphi  + Tpk*sphi)*g2             # ~zp
        # # + (Spk*sphi  - Tpk*cphi)*g2*k*z         # ~zp*z
        # # - (Sppk*cphi + Tppk*sphi)*k*zp          # ~zp*zp
        # # - (Sppk*sphi - Tppk*cphi)*k**2*z*zp     # ~z*zp*zp
        #   )
        # )
        # PoPi = om*qV0/(b**3*g*g2**2*gb3*m0c3)*((b**3*g*g2)*(
        #     Spk*cphi*g2 
        #     + Spk*g2*k*sphi*z 
        #     - Sppk*cphi*k*zp 
        #     # - Sppk*k**2*sphi*z*zp 
        #     - Tpk*cphi*g2*k*z 
        #     + Tpk*g2*sphi 
        #     # + Tppk*cphi*k**2*z*zp 
        #     - Tppk*k*sphi*zp)  
        # - 3*gb3*zp*(
        #     Spk*cphi*g2 
        #     # + Spk*g2*k*sphi*z 
        #     # - Sppk*cphi*k*zp 
        #     # - Sppk*k**2*sphi*z*zp 
        #     # - Tpk*cphi*g2*k*z 
        #     + Tpk*g2*sphi 
        #     # + Tppk*cphi*k**2*z*zp 
        #     # - Tppk*k*sphi*zp
        #     )) 
            # PoPi = om*qV0*(1/gb3 - 3*zp/(b**3*g*g2))*((Spk - Sppk*k*zp/g2)*(cphi + k*sphi*z) + (Tpk - Tppk*k*zp/g2)*(-cphi*k*z + sphi))/m0c3
        def OUTminusIN(qV0,Tk,Sk,Tpk,Spk,Tppk,Sppk,k,Db2b,Dphi,om,b,g,gb3,mc3,cphi,sphi):
            WoWi = qV0*(
            # Db2b*Dphi*Spk*cphi*k 
            # + Db2b*Dphi*Tpk*k*sphi 
            + Db2b*Spk*k*sphi 
            - Db2b*Tpk*cphi*k 
            - Dphi*Sk*cphi 
            - Dphi*Tk*sphi 
            - Sk*sphi 
            + Tk*cphi)
        
            PoPi = -om*qV0*(-3*Db2b*gb3*(
            # Dphi*Spk*sphi 
            # - Dphi*Tpk*cphi 
            - Spk*cphi 
            - Tpk*sphi) 
            + b**3*g*(Dphi*Spk*sphi 
            - Dphi*Tpk*cphi 
            - Spk*cphi 
            - Tpk*sphi))/(b**3*g*gb3*mc3)
            return (WoWi,PoPi)

        polies   = self.polies
        c        = PARAMS['clight']
        m0c2     = self.particle.m0c2
        m0c3     = m0c2*c
        omega    = self.omega
        ttf      = 0.
        matrix   = NP.eye(MDIM,MDIM)

        # initialise loop variables
        Ts       = self.particle.tkin                 # T-soll
        phis     = self.phisoll                       # phi-soll

        for poly in polies:   # each poly is a slice of the full Ez-dist
            # IN variables
            gammas_in  = 1. + Ts/m0c2               # gamma 
            gbs_in     = sqrt(gammas_in**2-1.)      # (gamma*beta)
            betas_in   = gbs_in/gammas_in           # beta

            ks     = omega/(c*betas_in)    # [1/m] omega/(beta*c)
            qV0    = self.V0(poly)         # [MV]
            Tks    = self.T(poly,ks)
            Sks    = self.S(poly,ks)
            Tpks   = self.Tp(poly,ks)
            Spks   = self.Sp(poly,ks)
            Tppks  = self.Tpp(poly,ks)  
            Sppks  = self.Spp(poly,ks) 

            sphis  = sin(phis)
            cphis  = cos(phis)

            # kin.energy and phase increase ref-particle
            cphi = cphis
            sphi = sphis
            om = self.omega
            bets = betas_in
            gams = gammas_in
            gb3s = (bets*gams)**3
            Db2b = 0.
            Dphi  = 0.
            (dTs,dPhis) = OUTminusIN(qV0,Tks,Sks,Tpks,Spks,Tppks,Sppks,ks,Db2b,Dphi,om,bets,gams,gb3s,m0c3,cphi,sphi)

            Ts_out    = Ts + dTs
            phis_out  = phis + dPhis
            ttf       = ttf + qV0

            # OUT variables
            gammas_out  = 1. + Ts_out/m0c2            # gamma 
            gbs_out     = sqrt(gammas_out**2-1.)      # (gamma*beta)
            #======================================================================================="""        
            # OXAL-matrix 
            (r44,r54) = OUTminusIN(qV0,Tks,Sks,Tpks,Spks,Tppks,Sppks,ks,betas_in,gammas_in,m0c3,phis,omega,z=1,zp=0)
            (r45,r55) = OUTminusIN(qV0,Tks,Sks,Tpks,Spks,Tppks,Sppks,ks,betas_in,gammas_in,m0c3,phis,omega,z=0,zp=1)
            #======================================================================================="""        
            # {z, dP/P}: linear sub matrix
            mx = NP.eye(MDIM,MDIM)
            mx[Ktp.z, Ktp.z] = r44 + 1; mx[Ktp.z, Ktp.zp ] = r45*(1/ks)
            mx[Ktp.zp,Ktp.z] = r54; mx[Ktp.zp, Ktp.zp] = r55*(1/ks) + 1
            # {x,x'}: linear sub-matrix
            factor2 = qV0*omega/(2.*m0c3*gbs_out*gbs_in**2)
            mx[Ktp.xp,Ktp.x ] = -factor2 * (Tks*sphis + Sks*cphis)
            mx[Ktp.xp,Ktp.xp] = gbs_in/gbs_out
            # {y,y'}: linear sub-matrix
            mx[Ktp.yp,Ktp.y]  = mx[Ktp.xp, Ktp.x]
            mx[Ktp.yp,Ktp.yp] = mx[Ktp.xp, Ktp.xp]
            # energy and length increase
            mx[Ktp.T,Ktp.dT] = dTs
            mx[Ktp.S,Ktp.dS] = 0     # 0 length: oxal-gap is kick

            # left multiplication of slice-matrix with oxal-matrix
            matrix = NP.matmul(mx,matrix)

            # refresh loop variables
            Ts    = Ts_out
            phis  = phis_out
        self.deltaW    = matrix[Ktp.T,Ktp.dT]
        self.ttf       = abs(self.deltaW/ttf)
        self.particlef = Proton(self.particle.tkin+self.deltaW)
        self.matrix    = matrix
        pass
        return
