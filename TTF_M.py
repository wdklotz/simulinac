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
from math import sin,cos,tan,pi,sqrt
import numpy as NP
import unittest
from setutil import PARAMS,I0,I1,MDIM,WConverter,Twiss,Proton,OutOfRadialBoundEx
from setutil import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,DSKOO
from setutil import DEBUG_ON,DEBUG_OFF,log_what_in_interval
from Ez0 import SFdata
from separatrix import w2phi

twopi = 2.*pi
pihalf = pi/2

def cot(x):
    return -tan(x+pihalf)

class TTF_G(IGap.IGap):
    """ 3 Point TTF Model (A.Shishlo/J.Holmes ORNL/TM-2015/247)"""
    def __init__(self):
        pass

    def configure(self,**kwargs):
        # static
        self.mapping      = 'ttf'
        self.kwargs       = kwargs
        self.label        = 'TTF' 
        # injected
        self.EzPeak    = kwargs.get('EzPeak',None)
        self.dWf       = kwargs.get('dWf',1)
        self.phisoll   = kwargs.get('phisoll',None)
        self.cavlen    = kwargs.get('cavlen',None)
        self.freq      = kwargs.get('freq',None)
        self.SFdata    = kwargs.get('SFdata',None)
        self.particle  = kwargs.get('particle',Proton(50.))
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)
        # calculated
        self.lamb      = PARAMS['clight']/self.freq
        self.gap       = self.cavlen*0.875   # 87.5%: of cavlen
        self.omega     = twopi*self.freq
        self.polies    = self.poly_slices()
        self.m0c2      = self.particle.m0c2
        self.m0c3      = self.m0c2**3
        self.ttf       = None
        self.matrix    = None
        self.deltaW    = None
        self.particlef = None
        self.master    = None
        pass

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return self.ttf_map(i_track)

    def toString(self):
        return f'{self.mapping} mapping in: TTF_M.ttf_map()'

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

    def adjust_energy(self,tkin):
        self.particle  = Proton(tkin)
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.deltaW    = self.EzPeak * self.ttf * self.gap * cos(self.phisoll)
        self.particlef = Proton(tkin + self.deltaW)
        self.T3D_matrix()
        pass

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
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz*cot(k*dz))
        t  = f1*f2
        DEBUG_OFF('TTF_G: (T,k) {}'.format((t,k)))
        return t
    def S(self, poly, k):    # A.Shishlo/J.Holmes (4.4.7)
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2       # [1/m] --> [1/cm]
        f1 = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.-k*dz*cot(k*dz)
        s  = f1*f2
        DEBUG_OFF('TTF_G: (S,k) {}'.format((s,k)))
        return s
    def Tp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.8)
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz*cot(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp
    def Sp(self, poly, k):   # A.Shishlo/J.Holmes (4.4.9)
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        sp  = 2*a*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        sp  = sp*(dz**2-2./k**2+dz*cot(k*dz)*2/k)
        sp  = sp*1.e-2     # [cm] --> [m]
        return sp
    def V0(self, poly):      # A.Shishlo/J.Holmes (4.4.3)
        E0 = poly.E0                          # [MV/m]
        b  = poly.b                           # [1/cm**2]
        dz = poly.dz                          # [cm]
        v0 = E0*(2*dz+2./3.*b*dz**3)*1.e-2    # [MV]
        return v0
    def poly_slices(self):
        """ Slice the RF gap """
        slices = []
        zl = -self.gap/2.*100.   # [m] --> [cm]
        zr = -zl
        for poly in self.SFdata.polies:
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slices.append(poly)
        return slices
    def ttf_map(self, i_track):
        def ttf_formeln(particle,phiIN,poly,r):
            omega= self.omega
            c    = PARAMS['clight']
            # m0c2 = self.m0c2
            m0c3 = self.m0c3
            g    = particle.gamma
            b    = particle.beta
            tkin = particle.tkin
            gb   = particle.gamma_beta
            gb3  = gb**3
            k    = self.omega/(c*b)
            V0   = self.V0(poly)
            Tk   = self.T(poly,k)
            Sk   = self.S(poly,k)
            Tkp  = self.Tp(poly,k)
            Skp  = self.Sp(poly,k)
            cphi = cos(phiIN)
            sphi = sin(phiIN)
            kr   = k/g*r
            i0   = I0(kr)
            i1   = I1(kr)
            DW   = V0*i0*(Tk*cphi - Sk*sphi)   # Shishlo 4.3.1
            Dphi = V0*omega/(m0c3*gb3)*(i0*(Tkp*sphi+Skp*cphi)+g*r*i1*(Tk*sphi+Sk*cphi)) # Shislo 4.3.2
            return (DW,Dphi,i0,i1,V0,Tk,Sk,Tkp,Skp)
        #===================== ttf_map =====================================
        # Teilchen Koordinaten IN
        x         = i_track[XKOO]       # [0]
        xp        = i_track[XPKOO]      # [1] Dx/Ds
        y         = i_track[YKOO]       # [2]
        yp        = i_track[YPKOO]      # [3] Dy/Ds
        z         = i_track[ZKOO]       # [4] z
        zp        = i_track[ZPKOO]      # [5] dp/p
        T         = i_track[EKOO]       # [6] kinetic energy Sollteilchen
        S         = i_track[SKOO]       # [8] position gap

        max_r     = 0.05                # max radial excursion [m]
        r         = sqrt(x**2+y**2)     # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(S)
            sys.exit()

        lamb         = self.lamb
        phisoll      = self.phisoll
        freq         = self.freq

        """
        ptcles Soll particle
        Ws     kin energy Soll
        phis   phase Soll

        ptcle Off particle
        W     kin energy 
        phi   phase
        """
        # Soll
        ptcles   = self.particle
        Ws       = ptcles.tkin
        phis     = self.phisoll

        # Off
        converter = WConverter(Ws,freq)
        Dphi      = converter.zToDphi(z)    # z->phi
        DW        = converter.Dp2pToDW(zp)  # zp->W
        phi       = phis + Dphi
        W         = Ws + DW
        ptcle     = Proton(W)

        self.deltaW  = 0.
        self.ttf     = 0.
        
        for poly in self.polies:
            """ Map through poly; use Formel 4.3.1 & 4.3.2 A.Shishlo/J.Holmes """
            gbIs = ptcles.gamma_beta
            cphi = cos(phis)
            sphi = sin(phis)

            # Soll IN
            (DWs,Dphis,i0,i1,V0,Tk,Sk,Tkp,Skp) = ttf_formeln(ptcles,phis,poly,r)
            # Soll OUT  (W,phi)
            Ws         = Ws + DWs
            phis       = phis + Dphis
            ptcles     = Proton(Ws)
            gbOs       = ptcles.gamma_beta

            # Off IN
            (DW,Dphi,i0x,i1x,V0x,Tkx,Skx,Tkpx,Skpx) = ttf_formeln(ptcle,phi,poly,r)
            # Off OUT
            W         = W + DW
            phi       = phi + Dphi
            ptcle     = Proton(W)

            i12r      = i1/r if r > 1.e-6 else pi/lamb/gbIs
            faktor    = V0/(self.m0c2*gbIs*gbOs)*i12r
            gbIs2gbOs = gbIs/gbOs
            xp        = gbIs2gbOs*xp-faktor*(Tk*sphi + Sk*cphi)*x
            yp        = gbIs2gbOs*yp-faktor*(Tk*sphi + Sk*cphi)*y

            self.deltaW   = self.deltaW + DWs
            self.ttf      = self.ttf + Tk
            pass
        
        T = T + self.deltaW
        # # Umwandlung longitudinale Koordinaten
        zO      = converter.DphiToz(phi-phis)
        zpO     = converter.DWToDp2p(W-Ws)
        f_track = NP.array([x,xp,y,yp,zO,zpO,T,1.,S,1.])
        self.particlef = ptcles
        self.ttf       = self.ttf/len(self.polies) # gap's ttf

        # S1 = 1    # from
        # S2 = 20    # to 
        # log_what_in_interval(S,(S1,S2),f'TTF_M.ttf_map: f_track: {f_track}\n')

        return f_track

def ttf(lamb, gap, beta):
    """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
    x = gap/(beta*lamb)
    res =NP.sinc(x)
    return res

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