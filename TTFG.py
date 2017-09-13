import matplotlib.pyplot as plt
import numpy as np
from math import sin,cos,tan,pi,exp,fabs,pow,sqrt,fmod,radians
from collections import namedtuple

from setutil import FLAGS,PARAMS,DEBUG,Proton
import elements as ELM
from Ez0 import SFdata,KpolySF,V0,Tk,Tkp,Sk,Skp

## Transition Time Factors RF Gap Model
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
        self.Ez0    = max([dpoint.Ez for dpoint in self.Ez])    # peak of EZ-values
        self.dWf    = dWf
        self.lamb   = PARAMS['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.anz    = 10                                        # nboff slices (imutable!)
        self.poly   = KpolySF(self.Ez,anz=self.anz)
        self.ks     = (2*pi)/(self.lamb*self.particle.beta)
        self.v0s()
    def v0s(self):
        zl          = -self.gap/2.*100.              #[cm]
        zr          = -zl     
        t     = Tk (self.poly,self.ks,zl,zr)
        V     = V0 (self.poly,zl,zr)                 #[MV] 
        self.V0 = [V[i]*t[i] for i in range(len(t))] # energy gain of sollteichen
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        DEBUG('TTFG: adapt_for_energy ',tkin)
        self.particle(tkin)
        self.v0s()
        return self
    def make_slices(self,anz=0):
        anz = self.anz      # prevent caller to set anz (shall stay fixed)
        slices = np.full(anz,self)
        self.active_slice = 0
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
    slices = ttfg.make_slices(anz=100)
    DEBUG('slices:\n',slices)
    DEBUG('ttfg.__dict__',ttfg.__dict__)
    
if __name__ == '__main__':
    test0()