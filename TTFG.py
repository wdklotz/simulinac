# import matplotlib.pyplot as plt
# import numpy as np
# from math import sin,cos,tan,pi,exp,fabs,pow,sqrt,fmod,radians
import sys
from math import sin,cos,tan,pi,radians,degrees,sqrt
from collections import namedtuple
from copy import copy

from setutil import FLAGS,PARAMS,DEBUG,Proton,I0,I1,tblprnt
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM
from Ez0 import SFdata

## Transition Time Factors RF Gap Model
class TTFGslice(object):
    def __init__(self,ttfg,polyval,particle):
        self.elemnt     = ttfg    # the element this slice is part off
        self.polyval    = polyval # ACHTUNG: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.particle   = copy(particle)   # incoming sync. particle
        self.beta       = self.particle.beta
        self.ks         = (2*pi)/(PARAMS['wellenlänge']*self.beta)*1.e-2    #[1/cm]
        self.Tk         = self._T (self.polyval,self.ks)
        self.Tkp        = self._Tp(self.polyval,self.ks)
        self.dWf        = self.elemnt.dWf
        self.V0         = self._V0(self.polyval)
        self.phi        = None
    def setPhase(self,value):
        self.phi = value
    def getPhase(self):
        return self.phi
    def getW(self):
        return self.particle.tkin
    def getDW(self):        # delta-Ws in this slice
        return self._DW()
    def getDPhase(self):     # delta-phi in this slice
        return self._DPhase()
    def getV0(self):
        return self.V0
    def getTk(self):
        return self.Tk
    def getTkp(self):
        return self.Tkp
    def getLen(self):
        return self.polyval.dz
    def _T(self,poly,k):
            a  = poly.a
            b  = poly.b
            dz = poly.dz
            f1 = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
            f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
            # DEBUG('Tk: (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
            t  = f1*f2
            return t
    def _Tp(self,poly,k):
            a   = poly.a
            b   = poly.b
            dz  = poly.dz
            tp  = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
            tp  = tp * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
            return tp
    def _V0(self,poly):
            E0 = poly.E0                          #[MV/m]
            b  = poly.b
            dz = poly.dz                          #[cm]
            v0 = E0*(2*dz+2./3.*b*dz**3)*1.e-2    #[MV]
            v0 = v0*self.dWf
            return v0
    def _DW(self):
        V0    = self.V0
        Tk    = self.Tk
        phase = self.phi
        res   = V0*Tk*cos(phase)
        return res
    def _DPhase(self):
        f     = PARAMS['frequenz']
        c     = PARAMS['lichtgeschwindigkeit']
        m0c2  = self.particle.e0
        bg    = self.particle.gamma_beta
        phase = self.phi
        # pvlen = self.polyval.zl-self.polyval.zr  # pvlen = z(t)-z(t+dt) < 0
        pvlen = -self.polyval.dz
        res   = self.V0*2*pi*f/(m0c2*c)/bg**3
        res   = res * self.Tkp*sin(phase)
        res   = res*pvlen*1.e-2    #[cm]-->[m]
        return res
    def adapt_for_energy(self,tkin):
            self.particle(tkin)
            self.beta    = self.particle.beta
            self.ks      = (2*pi)/(PARAMS['wellenlänge']*self.beta)*1.e-2    #[1/cm]
            self.Tk      = self._T (self.polyval,self.ks)
            self.Tkp     = self._Tp(self.polyval,self.ks)
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
        self.dWf    = dWf

        self.lamb   = PARAMS['wellenlänge']
        self.slices = self._make_slices()
        self.deltaW = None       # delta-W of synchronous particle through  TTF-gap

    def shorten(self,l=0.):
        # no shortening for TTFG yet
        return self
    def adapt_for_energy(self,tkin):
        next_phase = self.phis
        next_tkin  = tkin
        for cnt,slice in enumerate(self.slices):
            DEBUG('TTFGslice[{}]: phase {},\ttkin {}'.format(cnt,degrees(next_phase),next_tkin))
            slice.setPhase(next_phase)
            slice.adapt_for_energy(next_tkin)
            DPhase     = slice.getDPhase()
            DW         = slice.getDW()
            next_phase = slice.getPhase()+DPhase
            next_tkin  = slice.particle.tkin+DW
        self.deltaW = self.slices[-1].particle.tkin-self.slices[0].particle.tkin
        return self
    def _make_slices(self,anz=0):
        slices = []
        zl = -self.gap/2.*100.   # [m] --> [cm]
        zr = -zl
        for poly in self.Ez.Ez_poly():
            zil = poly.zl
            zir = poly.zr
            if zil < zl or zir > zr: continue
            slice = TTFGslice(self,poly,self.particle)  # instanciate TTFGslices
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
        Ti        = i_track[EKOO]       # [6] summe aller deltaW
        # DTi       = i_track[DEKOO]    # [7] 1 
        si        = i_track[SKOO]       # [8] summe aller laengen
        # Dsi       = i_track[LKOO]     # [9] 1 

        m0c2       = self.particle.e0
        lamb       = self.lamb
        twopi      = 2.*pi
        c          = PARAMS['lichtgeschwindigkeit']

        dbTab1Rows  = []
        dbTab2Rows  = []
        
        for slice in self.slices:
            dbTab1Row = []
            particlei = slice.particle
            betai     = particlei.beta
            gammai    = particlei.gamma
            gbi       = particlei.gamma_beta
            qV0       = slice.getV0()
            Tki       = slice.getTk()
            Tkpi      = slice.getTkp()
            Wsi       = slice.getW()         # sync. part. W @ (in)
            phisi     = slice.getPhase()     # sync. part. phi @ (in)
            Ski       = 0.
            Skpi      = 0.

            r = sqrt(xi**2+yi**2)            # radial coordinate
            Kr = (twopi*r)/(lamb*gbi)
            i0 = I0(Kr)                      # bessel function
            i1 = I1(Kr)                      # bessel function

            conv_zDPhi  = -twopi/(betai*lamb)       # conversion factor z --> Dphi @ (in)
            conv_dPdT   = 1./(m0c2*betai**2*gammai) # conversion factor dp/p --> dT @ (in)

            Dphii   = conv_zDPhi*zi          # zi -> dPhi
            phii    = phisi+Dphii            # part. phi @ (in)
            commonf = twopi*qV0/(lamb*m0c2*c*gbi**3)
            phifmphii = commonf*(i0*r*(Tkpi*sin(phii)+Skpi*cos(phii))+
                gammai*r*i1*(Tki*sin(phii)+Ski*cos(phii)))*zi # Formel 4.3.2 A.Shishlo
            phisf = phisi+slice.getDPhase()
            phisfmphisi = phisf - phisi
            Dphi = phifmphii - phisfmphisi

            DWi     = conv_dPdT*zpi          # dp/p --> dT @ (in)
            Wi      = Wsi + DWi              # part. W @ (in)
            if Wi < 0.:
                raise RuntimeError('negative kinetic energy {:8.4g}'.format(Wi))
                sys.exit(1)
            Wsf = Wsi+slice.getDW()          # sync. part W @ (f)
            particlef = copy(particlei)(tkin=Wsf)
            betaf = particlef.beta
            gammaf = particlef.gamma
            gbf = particlef.gamma_beta

            conv_Dphiz = -(betaf*lamb)/twopi    # conversion factor Dphi --> z @ (f)
            conv_dTdP  = (m0c2*betaf**2*gammaf) # conversion factor dT --> dp/p @ (f)
            WsfmWsi = Wsf - Wsi
            WfmWi= qV0*i0*(Tki*cos(phii)-Ski*sin(phii))  # Formel 4.3.1 A.Shishlo     
            DWf = WfmWi - WsfmWsi
            
            zf  = conv_Dphiz*Dphi            # Dphf --> zf
                                             # zf   = gbf/gbi*zi  ?? needed ??
            zpf = conv_dTdP*DWf              # dT --> dp/p

            xf = xi                          # no change
            yf = yi                          # no change

            commonf = qV0/(m0c2*gbi*gbf)*i1                   # common factor
            xpf = gbi/gbf*xpi-xi/r*commonf*(Tki*sin(phii)+Ski*cos(phii)) # Formel 4.3.3 A.Shishlo
            ypf = gbi/gbf*ypi-yi/r*commonf*(Tki*sin(phii)+Ski*cos(phii)) # Formel 4.3.3 A.Shishlo

            Tf   = Ti + slice.getDW()
            sf   = si + slice.getLen()*1.e-2                  # [cm] --> [m]

            dbTab1Row = [
                '{:8.4f}'.format(Wsi),
                '{:8.4f}'.format(r),
                '{:8.4f}'.format(Kr),
                '{:8.4f}'.format(i0),
                '{:8.4f}'.format(i1),
                '{:8.4g}'.format(degrees(Dphii))
                ]
            dbTab1Rows.append(dbTab1Row)
            dbTab2Row = [
                '{:8.4g}'.format(zf*1.e3),
                '{:8.4g}'.format(degrees(Dphi)),
                '{:8.4g}'.format(zpf*1.e2),
                '{:8.4g}'.format(DWf*1.e3),
                '{:8.4g}'.format(xpf*1.e3),
                '{:8.4g}'.format(ypf*1.e3),
                ]
            dbTab2Rows.append(dbTab2Row)
#   f_track = NP.array([xf,xpf,yf,ypf,zf,zpf,DWs,DTf,sf,Dsf])
        dbTab1Headr = ['Ws(in)','r','Kr','i0','i1','Dphi(in)']
        dbTab2Headr = ['z[mm]','Dphi[deg]','dp/p[%]','DWf[KeV]','xp[mrad]','yp[mrad]']
        DEBUG(' at position IN:\n'+(tblprnt(dbTab1Headr,dbTab1Rows)))
        DEBUG(' at position OUT:\n'+(tblprnt(dbTab2Headr,dbTab2Rows)))
        pass
def test0():
    from tracks import Track
    import numpy as np
    
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    
    ttfg = TTFG(length=0.044,gap=0.048,Ez=SF_tab)
    ttfg.adapt_for_energy(50.)
    DEBUG('ttfg.__dict__',ttfg.__dict__)
    slices = ttfg.slices
    for slice in slices:
        DEBUG('slice\n',slice.__dict__)

    z = 1.e-3
    x=y=1.e-2
    tkin=10.
    # a track-point is that:       x   x'  y  y'  z   z'  W    1   s   1
    track = Track(start=np.array([ x,  0., y, 0., z,  0.,tkin, 1., 0., 1.]))
    # track point @ (i)
    track_i = track.first()
    ttfg.map(track_i)

if __name__ == '__main__':
    test0()