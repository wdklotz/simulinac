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

#TODO: Formeln für Tpp und Spp
#TODO: unittests
"""
import sys
import IGap
import unittest
import numpy as NP
from math import sin,cos,tan,sqrt,pi,degrees,radians
from setutil import FLAGS,PARAMS,Ktp,MDIM,Proton,DEBUG_ON,DEBUG_OFF,Proton
from setutil import wrapRED,mxprnt,Twiss,WConverter,dictprnt
from Ez0 import SFdata
from separatrix import w2phi
import Ez0 as EZ

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
        pass

    def configure(self,**kwargs):
        self.dWf          = FLAGS['dWf']
        self.mapping      = 'oxal'
        self.kwargs       = kwargs
        self.label        = 'OX'

        self.EzPeak    = kwargs.get('EzPeak',None)
        self.phisoll   = kwargs.get('phisoll',None)
        self.cavlen    = kwargs.get('cavlen',None)
        self.freq      = kwargs.get('freq',None)
        self.SFdata    = kwargs.get('SFdata',None)
        self.particle  = kwargs.get('particle',Proton(50.))
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)
        self.sec       = kwargs.get('sec')

        self.omega     = twopi*self.freq
        self.lamb      = PARAMS['clight']/self.freq
        self.ttf       = None
        self.deltaW    = None
        self.particlef = None
        self.matrix    = None
        self.master    = None
        self.polies    = self.poly_slices()
        self.OXAL_matrix_1()
        pass

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return NP.dot(self.matrix,i_track)

    def toString(self):
        return mxprnt(self.matrix,'4g')

    def isAccelerating(self):
        return True

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
        self.particle = Proton(tkin)
        self.OXAL_matrix()
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
        self.matrix    = matrix
        self.particlef = Proton(self.particle.tkin+self.deltaW)
        return
    def OXAL_matrix_1(self):
        #====================================================================== sympy ==========
        def OUTminusIN(qV0,T,S,Tp,Sp,Tpp,Spp,k,b,g,m0c3,phis,om,z=0,zp=0):
            """ k = omega/(c*beta), b = beta, g = gamma, om=omega
                T,S,Tp,Sp,Tpp,Spp = T(k),S(k),T'(k),S'(k),T''(k),S''(k) Fourierfaktoren
            """
            g2   = g**2      # gamma**2
            cphi = cos(phis)
            sphi = sin(phis)  
            Db2b = zp/g2     # Delta-beta/beta
            Dp2p = - k*z     # Delta-phi/phi

            WoWi = qV0/g2*(
            #  - (Spks*cos(p) + Tpks*k*sin(p))*Dp2p*k**2*z    # O2
               + (Sp*sphi - Tp*cphi)*Dp2p*k           # O1   delta-p/p
               + (S*cphi  + T*sphi)*g2*k*z          # O1   z
               - (S*sphi  - T*cphi)*g2)             # O0

            PoPi = -om*qV0/(b**3*g**5*m0c3)*(3*Dp2p*(
            #	- (Sppks*sin(p) - Tppks*cos(p))*Dp2p*k**2*z    # O3
            #	- (Sppks*cos(p) + Tppks*sin(p))*Dp2p*k         # O2
            #	+ (Spks*sin(p)  - Tpks*cos(p))*g**2*k*z        # O2
                + (Sp*cphi  + Tp*sphi)*g2)           # ~Dp2p
            - 1*(
            #	- (Sppks*sin(p) - Tppks*cos(p))*Dp2p*k**2*z    # O2
                - (Spp*cphi + Tpp*sphi)*Dp2p*k         # ~Dp2p
                + (Sp*sphi  - Tp*cphi)*g2*k*z        # ~z
                + (Sp*cphi  + Tp*sphi)*g2)           # O0
            )
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
            (dTs,dPhis) = OUTminusIN(qV0,Tks,Sks,Tpks,Spks,Tppks,Sppks,ks,betas_in,gammas_in,m0c3,phis,omega,z=0,zp=0)

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
            mx[Ktp.z, Ktp.z] = r44; mx[Ktp.z, Ktp.zp ] = r45 
            mx[Ktp.zp,Ktp.z] = r54; mx[Ktp.zp, Ktp.zp] = r55 
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
        self.matrix    = matrix
        self.particlef = Proton(self.particle.tkin+self.deltaW)
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
        print(wrapRED('------------------ test_OXAL'))
        gap_parameter = dict(
            EzPeak    = 1,
            phisoll   = radians(-30.),
            gap       = 0.022,
            cavlen    = 0.44,
            freq      = 750e6,
        )
        fieldtab  = 'SF/CAV-FLAT-R135-L31.TBL'
        EzPeak    = gap_parameter['EzPeak']
        cavlen    = gap_parameter['cavlen']
        sfdata    = EZ.SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,L=cavlen/2.*100.)   # scaled field distribution
        gap_parameter['SFdata'] = sfdata

        oxag = OXAL_G()    # create object instance
        oxag.configure(**gap_parameter)
        print(f'transfer matrix OXAL for {oxag.particle.tkin} MeV (default)')
        print(oxag.toString())

        tkin = 100.
        print(f'transfer matrix OXAL for {tkin} MeV')
        oxag.adjust_energy(tkin)
        print(oxag.toString())
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
        res = oxag.waccept()
        dictprnt(what=res,text='waccept',njust=25)

        # injection_energies = [6,12,24,50,100,150,200]
        # EzPeak = 2.0; phisoll = radians(-30.); gap = 0.044; freq = 750.e6; fieldtab="SF/SF_WDK2g44.TBL"
        # gap_cm = 100.*gap   # gap in cm
        # sfdata = SFdata.InstanciateAndScale(fieldtab,EzPeak=EzPeak,IntgInterval=gap_cm)
        # EzAvg = sfdata.EzAvg
        # ID = "OXAL"
        # for injection_energy in injection_energies:
        #     oxal = OXAL_G(ID,EzAvg,phisoll,gap,freq,SFdata=sfdata,particle=Proton(injection_energy))
        #     tvectori = NP.array([0, 0, 0, 0, 0, 0, injection_energy, 1, 0, 1])
        #     tvectoro = NP.dot(oxal.matrix,tvectori)
        #     print(f'EzPeak={EzPeak}, gap={gap}, freq={freq*1e-6}, Wi={tvectori[6]:3.3f}, Wo={tvectoro[6]:3.3f}, dW={oxal.matrix[6,7]:.3e}')
        # return

if __name__ == '__main__':
    unittest.main()
