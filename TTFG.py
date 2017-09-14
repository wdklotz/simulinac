# import matplotlib.pyplot as plt
# import numpy as np
# from math import sin,cos,tan,pi,exp,fabs,pow,sqrt,fmod,radians
from math import sin,cos,tan,pi,radians,degrees
from collections import namedtuple
from copy import copy

from setutil import FLAGS,PARAMS,DEBUG,Proton
import elements as ELM
from Ez0 import SFdata

## Transition Time Factors RF Gap Model
class TTFGslice(object):
    def __init__(self,ttfg,poly,particle):
        self.driver     = ttfg
        self.polyval    = poly  #ACHTUNG: Polynomkoeffizienten in [1/cm] u. [1/cm**2]
        self.particle   = copy(particle)           # incoming sync. particle
        self.beta       = self.particle.beta
        self.ks         = (2*pi)/(PARAMS['wellenlänge']*self.beta)*1.e-2    #[1/cm]
        self.Tk         = self._T (self.polyval,self.ks)
        self.Tkp        = self._Tp(self.polyval,self.ks)
        self.V0         = self._V0(self.polyval)
        self.phi      = 0.      # ? radians(-25.) oder 0. ?
        # self.phi        = self.driver.phis
        self.DWs        = self._DWs()
        self.Dphase     = self._DPhase()
    def setPhase(self,value):
        self.phi = value
    def getPhase(self):
        return self.phi
    def getDWs(self):   # delta-Ws
        return self.DWs
    def getDphase(self): # delta-phi
        return self.Dphase
    def _T(self,poly,k):
            a  = poly.a
            b  = poly.b
            dz = poly.dz
            f1 = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
            f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
            # DEBUG('Tk: (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
            t = f1*f2
            return t
    def _Tp(self,poly,k):
            a   = poly.a
            b   = poly.b
            dz  = poly.dz
            tp = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
            tp = tp * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
            return tp
    def _V0(self,poly):
            E0 = poly.E0                          #[MV/m]
            b  = poly.b
            dz = poly.dz                          #[cm]
            v0  = E0*(2*dz+2./3.*b*dz**3)*1.e-2   #[MV] 
            return v0
    def _DWs(self):
        V0    = self.V0
        Tk    = self.Tk
        phase = self.phi
        res = V0*Tk*cos(phase)
        return res
    def _DPhase(self):
        f     = PARAMS['frequenz']
        c     = PARAMS['lichtgeschwindigkeit']
        m0c2  = self.particle.e0
        bg    = self.particle.gamma_beta
        phase = self.phi
        res = self.V0*2*pi*f/(m0c2*c)/bg**3
        res = res * self.Tkp*sin(phase)
        return res
    def adapt_for_energy(self,tkin):
            self.particle(tkin)
            self.beta    = self.particle.beta
            self.ks      = (2*pi)/(PARAMS['wellenlänge']*self.beta)*1.e-2    #[1/cm]
            self.Tk      = self._T (self.polyval,self.ks)
            self.Tkp     = self._Tp(self.polyval,self.ks)
            self.DWs     = self._DWs()
            self.Dphase  = self._DPhase()
            return self
class TTFG(ELM.D):
    """
    Transition Time Factors RF Gap Model (A.Shishlo ORNL/TM-2015/247)
    """
    def __init__(self,
                    PhiSoll    = radians(PARAMS['soll_phase']),
                    fRF        = PARAMS['frequenz'],
                    label      = 'TTFG',
                    particle   = PARAMS['sollteilchen'],
                    gap        = PARAMS['spalt_laenge'],
                    length     = 0.,
                    Ez         = None,                # return of SFdata
                    dWf        = FLAGS['dWf']):
        super().__init__(label=label, particle=particle)
        if Ez == None: raise RuntimeError('TTFG: missing E(z) table')
        self.viseo  = 0.25
        self.phis   = PhiSoll    # [radians] soll phase
        self.freq   = fRF        # [Hz]  RF frequenz
        self.label  = label
        self.gap    = gap        # [m]
        self.length = length     # [m]
        self.Ez     = Ez         # the SF-table, a list(DataPoints)
        # self.Ez0    = max([dpoint.Ez for dpoint in self.Ez.Ez_table()])    # peak of EZ-values
        self.dWf    = dWf
        self.lamb   = PARAMS['wellenlänge']
        self.slices = self._make_slices()
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        DEBUG('TTFG: adapt_for_energy ',tkin)
        slicep = self.slices[0]                    # 1st slide
        dz = self.Ez.Ez_poly()[0].dz*1.e-2         # [cm] --> [m]
        Dphase = -2*pi/(self.particle.beta*self.lamb)*dz   # Achtung: (+) oder (-)?
        next_phase = slicep.getPhase()+Dphase
        slicep.setPhase(next_phase)
        slicep.adapt_for_energy(tkin)
        DEBUG('next_phase {}, next_tkin {}'.format(degrees(next_phase),tkin))
        Dphase = slicep.getDphase()
        DWs    = slicep.getDWs()    # 1st slide
        for slice in self.slices[1:]:
            next_phase = slicep.getPhase()+Dphase
            slice.setPhase(next_phase)             # increment phase (must be before energy)
            next_tkin = slicep.particle.tkin+DWs
            DEBUG('next_phase {}, next_tkin {}'.format(degrees(next_phase),next_tkin))
            slice.adapt_for_energy(next_tkin)      # increment energy
            Dphase = slice.getDphase()
            DWs    = slice.getDWs()
            slicep = slice
        return self
    def _make_slices(self,anz=0):
        slices = []
        zl = -self.gap/2.*100.   # [m] --> [cm]
        zr = -zl
        zl2zr = (zl,zr)
        poly = self.Ez.Ez_poly()
        for i in range(len(poly)):
            zil = poly[i].zl
            zir = poly[i].zr
            if zil < zl or zir > zr: continue
            # DEBUG('poly[{}]'.format(i))
            slice = TTFGslice(self,poly[i],self.particle)
            slices.append(slice)
        return slices
    def map(self,i_track):
        """
        Mapping of track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
        """
        xi        = i_track[XKOO]       # [0]
        xpi       = i_track[XPKOO]      # [1]
        yi        = i_track[YKOO]       # [2]
        ypi       = i_track[YPKOO]      # [3]
        zi        = i_track[ZKOO]       # [4] z-z0
        zpi       = i_track[ZPKOO]      # [5] dp/p - (dp/p)0
        Ti        = i_track[EKOO]       # [6] summe aller delta-T
        DTi       = i_track[DEKOO]      # [7] 1 immer
        si        = i_track[SKOO]       # [8] summe aller laengen
        Dsi       = i_track[LKOO]       # [9] 1 immer

        T          = self.tr
        # qE0L       = self.u0
        # qE0LT      = qE0L*T
        m0c2       = self.particle.e0
        lamb       = self.lamb
        phis       = self.phis
        twopi      = 2.*pi

       ##   particlesi = self.particle
        # tkinsi     = particlesi.tkin
        # betasi     = particlesi.beta
        # gammasi    = particlesi.gamma
        # gbsi       = particlesi.gamma_beta
        # 
        # DWs        = qE0LT*cos(phis)                # energy increase sync. particle
        # Wsf        = tkinsi + DWs 
        # particlesf = copy(particlesi)(tkin=Wsf)
        # betasf     = particlesf.beta
        # gammasf    = particlesf.gamma
        # gbsf       = particlesf.gamma_beta
        # 
        # condPdT = m0c2*betasi**2*gammasi
        # condTdP = 1./(m0c2*betasf**2*gammasf)

       ##   DWi       = condPdT*zpi                                # dp/p --> dT
        # Wi        = tkinsi + DWi    
        # if Wi < 0.:
        #     raise RuntimeError('negative kinetic energy {:8.4g}'.format(Wi))
        #     sys.exit(1)
        # particlei = copy(particlesi)(tkin=Wi)
        # betai     = particlei.beta
        # gbi       = particlei.gamma_beta
        # 
        # r = sqrt(xi**2+yi**2)                                 # radial coordinate
        # Kr = (twopi*r)/(lamb*gbi)
        # # DEBUG('r {:8.4g} Kr {:8.4g} gbi {:8.4g} Wi {:8.4g}'.format(r,Kr,gbi,Wi))
        # i0 = I0(Kr)                                           # bessel function
        # i1 = I1(Kr)                                           # bessel function

       ##   # THE MAP
        # zf      = gbsf/gbsi*zi
        # Dphii   = -zi*twopi/(betai*lamb)                       # z -> dPhi
        # phii    = Dphii+phis                                   # phi(in)
        # WfmWi   = qE0LT*i0*cos(phii)                           # W(f) - W(i)
        # WsfmWsi = DWs                                          # Ws(f) - Ws(i) particle energy increase
        # DWf     = WfmWi - WsfmWsi + DWi                        # DW(f)
        # zfp     = DWf*condTdP                                  # dT --> dp/p
        # 
        # xf   = xi     # x does not change
        # yf   = yi     # y does not change
        # DTf  = DTi    # 1
        # sf   = si     # because self.length always 0
        # Dsf  = Dsi    # 1
        # commonf = qE0LT/(m0c2*gbsi*gbsf)*i1                   # common factor
        # xpf  = gbsi/gbsf*xpi - xi/r*commonf*sin(phii)         # tranverse coordinate
        # ypf  = gbsi/gbsf*ypi - yi/r*commonf*sin(phii)         # tranverse coordinate

       ##   f_track = NP.array([xf,xpf,yf,ypf,zf,zfp,DWs,DTf,sf,Dsf])
        # # Wf      = Wi+DWf
        # # Wsi     = tkinsi
        # # DEBUG('RFB.rfb_map',
        # #     dict(
        # #         zi=zi, zpi=zpi, zf=zf, zfp=zfp,
        # #         Wi=Wi,Wf=Wf,Wsi=Wsi,Wsf=Wsf,phis=degrees(phis),phii=degrees(phii)
        # #         )
        # #     )  
        # return f_track
        pass
def test0():
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    Ez0_tab = SFdata(input_file,Epeak)
    
    ttfg = TTFG(length=0.044,gap=0.048,Ez=Ez0_tab)
    # DEBUG('ttfg.__dict__',ttfg.__dict__)
    ttfg.adapt_for_energy(100.)
    slices = ttfg.slices
    for slice in slices:
        DEBUG('slice\n',slice.__dict__)
    DEBUG('ttfg.__dict__',ttfg.__dict__)
    
if __name__ == '__main__':
    test0()