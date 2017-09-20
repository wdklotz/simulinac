# import matplotlib.pyplot as plt
# import numpy as np
# from math import sin,cos,tan,pi,exp,fabs,pow,sqrt,fmod,radians
import sys
from math import sin,cos,tan,radians,degrees,sqrt
from math import pi as PI
from collections import namedtuple
from copy import copy
import numpy as np

from setutil import FLAGS,PARAMS,DEBUG,Proton,I0,I1,tblprnt
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM
from Ez0 import SFdata

## DEBUG MODULE
def DEBUG_ON(string,arg='',end='\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg='',end='\n'):
    pass
DEBUG_MODULE   = DEBUG_OFF
DEBUG_MAP      = DEBUG_OFF
DEBUG_SOLL_MAP = DEBUG_OFF
twopi          = 2*PI

## Transition Time Factors RF Gap Model
class TTFGslice(object):
    def __init__(self,ttfg,polyval,particle):
        self.elemnt     = ttfg    # the element this slice is part off
        self.polyval    = polyval # polynom interval: ACHTUNG: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.particle   = copy(particle)   # incoming sync. particle
        self.beta       = self.particle.beta
        self.ks         = twopi/(self.elemnt.lamb*self.beta)*1.e-2    #[1/m] -->[1/cm]
        self.Tk         = self._T (self.polyval,self.ks)
        self.Tkp        = self._Tp(self.polyval,self.ks)
        self.dWf        = self.elemnt.dWf
        self.V0         = self._V0(self.polyval)
    def setPhase(self,value):
        self.phi = value   # phase of sync.particle
    def getPhase(self):
        try:
            return self.phi
        except:
            raise RuntimeError("phase @ slice not known yet!")
            sys.exit(1)
    def setDW(self,value):
        self.deltaw = value
    def getW(self):
        return self.particle.tkin
    def getDW(self):         # my deltaW
        return self._DW()
    def getDPhase(self):     # my delta-phi
        return self._DPhase()
    def getV0(self):
        return self.V0
    def getTk(self):
        return self.Tk
    def getTkp(self):
        return self.Tkp
    def _T(self,poly,k):
            a  = poly.a
            b  = poly.b
            dz = poly.dz
            f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
            f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
            t  = f1*f2
            return t
    def _Tp(self,poly,k):
            a   = poly.a
            b   = poly.b
            dz  = poly.dz
            tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
            tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
            tp  = tp*1.e-2     # [cm] --> [m]
            return tp
    def _V0(self,poly):
            E0 = poly.E0                          #[MV/m]
            b  = poly.b                           #[1/cm**2]
            dz = poly.dz                          #[cm]
            v0 = (2*dz+2./3.*b*dz**3)*1.e-2       #[m]
            v0 = v0*E0*self.dWf
            return v0
    def _DW(self):
        V0    = self.V0
        Tk    = self.Tk
        phase = self.getPhase()
        res   = V0*Tk*cos(phase)
        return res
    def _DPhase(self):
        pvlen = -self.polyval.dz
        res = self.ks*pvlen
        return res
    def adapt_for_energy(self,tkin):
            self.particle(tkin)
            self.beta    = self.particle.beta
            self.ks      = twopi/(self.elemnt.lamb*self.beta)*1.e-2    #[1/cm]
            self.Tk      = self._T (self.polyval,self.ks)
            self.Tkp     = self._Tp(self.polyval,self.ks)
            return self
class TTFG(ELM.I):
    """
    Transition Time Factors RF Gap Model (A.Shishlo ORNL/TM-2015/247)
    """
    def __init__(self,
                    PhiSoll    = radians(PARAMS['soll_phase']),
                    fRF        = PARAMS['frequenz'],
                    label      = 'TTFG',
                    particle   = PARAMS['sollteilchen'],
                    gap        = PARAMS['spalt_laenge'],
                    Ez         = None,                # return of SFdata
                    dWf        = FLAGS['dWf']):
        super().__init__(label=label,viseo=0.25,particle=particle)
        if Ez == None: raise RuntimeError('TTFG: missing E(z) table')
        self.phis   = PhiSoll    # [radians] soll phase
        self.freq   = fRF        # [Hz]  RF frequenz
        self.gap    = gap        # [m]
        self.dWf    = dWf
        self.lamb   = PARAMS['wellenlänge']
        self.Ez     = Ez         # the SF-table, a list-of-DataPoints
        try:
            self.Epeak  = Ez.Epeak
            self.slices = self._make_slices()
        except:
            raise RuntimeError("can't create object without SF-data!")
            sys.exit(1)
    def _make_slices(self):
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
    def adapt_for_energy(self,tkin):
        def adapt_slices_for_energy(tkin):
            next_phase = self.phis
            next_tkin  = tkin
            Tklist = []
            for slice in self.slices:
                slice.setPhase(next_phase)
                slice.adapt_for_energy(next_tkin)
                DPhase     = slice.getDPhase()
                DW         = slice.getDW()
                Tklist.append(slice.getTk())
                next_phase = slice.getPhase()+DPhase
                next_tkin  = slice.particle.tkin+DW
            deltaW   = next_tkin-tkin
            self.tr  = min(Tklist)
            return deltaW
        # body --------------------------------------------------------------        
        self.__init__(PhiSoll     = self.phis,
                      fRF         = self.freq,
                      label       = self.label,
                      particle    = self.particle(tkin),
                      gap         = self.gap,
                      Ez          = self.Ez,
                      dWf         = self.dWf)

        self.deltaW = adapt_slices_for_energy(tkin)  # deltaW for lin.map
        self.matrix[EKOO,DEKOO] = self.deltaW
        return(self)
    def make_slices(self,anz=0):
        return [self]
    def shorten(self,l=0.):
        return self
    def map(self,i_track):
        def slice_map(track):
            """
            Mapping of track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
            """
            def Dphiout_phiin(conv,g,r,i0,i1,tkp,skp,tk,sk,phi):
                # Formel 4.3.2 A.Shishlo
                return  conv*i0*(tkp*sin(phi)+skp*cos(phi)+g*r*i1*(tk*sin(phi)+sk*cos(phi)))
            def Dwout_win(conv,i0,tk,sk,phi):
                # Formel 4.3.1 A.Shishlo
                return conv*i0*(tk*cos(phi)-sk*sin(phi))
            def z2phase(conv,z,ps):
                return conv*z+ps
            def phase2z(conv,p,ps):
                return (p-ps)*conv
            def zp2W(conv,zp):
                return conv*zp
            def W2zp(conv,W,Ws):
                return conv*(W-Ws)

            x        = track[XKOO]       # [0]
            xp       = track[XPKOO]      # [1]
            y        = track[YKOO]       # [2]
            yp       = track[YPKOO]      # [3]
            z        = track[ZKOO]       # [4] z~(phi-phis)
            zp       = track[ZPKOO]      # [5] dp/p~dT
            T        = track[EKOO]       # [6] summe aller dT
            s        = track[SKOO]       # [8] summe aller laengen

            DEBUG_MAP("\n(x,x',y,y') ",arg=(x,xp,y,yp))
            DEBUG_MAP("(z,z') ",arg=(z,zp))

            m0c2       = self.particle.e0
            lamb       = self.lamb
            c          = PARAMS['lichtgeschwindigkeit']
   
            dbTab1Rows  = []      # for DEBUGGING
            dbTab2Rows  = []      # for DEBUGGING

            particlei = slice.particle
            betai     = particlei.beta
            gammai    = particlei.gamma
            gbi       = particlei.gamma_beta
            qV0       = slice.getV0()
            Tk        = slice.getTk()
            Tkp       = slice.getTkp()
            pINs      = slice.getPhase()     # soll phase  @ (i)
            WINs      = slice.getW()         # soll energy @ (i)

            r = sqrt(x**2+y**2)              # radial coordinate
            Kr = (twopi*r)/(lamb*gbi)
            i0 = I0(Kr)                      # bessel function
            i1 = I1(Kr)                      # bessel function

            DEBUG_MAP('(r,Kr,i 0,i1)',(r,Kr,i0,i1))

            conv_z2phase = -twopi/(betai*lamb)       # conversion-factor z --> (phi-phis)
            pIN = z2phase(conv_z2phase,z,pINs)       # phase(IN)
            common = twopi*qV0/(lamb*m0c2*gbi**3)
            dps = Dphiout_phiin(common,gammai,0.,1.,0.,Tkp,0.,Tk,0.,pINs)#pf-pi soll
            dp  = Dphiout_phiin(common,gammai,r ,i0,i1,Tkp,0.,Tk,0.,pIN) #pf-pi teilchen
            pOUTs = pINs+dps
            pOUT  = pIN+dp

            DEBUG_MAP('(z,pINs,pIN,dps,dp,pOUTs,pOUT)',(z,pINs,pIN,dps,dp,pOUTs,pOUT))

            conv_zp2W = (m0c2*betai**2*gammai)       # conversion-factor dp/p --> dW
            WIN = zp2W(conv_zp2W,zp)                 # energy(IN)
            dws = Dwout_win(qV0,1.,Tk,0.,pINs)
            dw  = Dwout_win(qV0,i0,Tk,0.,pIN)
            WOUTs = WINs+dws
            WOUT  = WIN+dw
    
            particlef = copy(particlei)(tkin=WOUTs)
            betaf = particlef.beta
            gammaf = particlef.gamma
            gbf = particlef.gamma_beta

            DEBUG_MAP('(WINs,WOUTs,WIN,WOUT,dws,dw)',(WINs,WOUTs,WIN,WOUT,dws,dw))

            common = qV0/(m0c2*gbi*gbf)*i1
            if r > 0.:
                xpf = gbi/gbf*xp-x/r*common*(Tk*sin(pIN)) # Formel 4.3.3 A.Shishlo (Sk=0.)
                ypf = gbi/gbf*yp-y/r*common*(Tk*sin(pIN))
            elif r == 0.:
                xpf = gbi/gbf*xp
                ypf = gbi/gbf*yp

            conv_phase2z = -(betaf*lamb)/twopi       # conversion-factor phi-phis --> z 
            conv_W2zp    = 1./(m0c2*betaf**2*gammaf) # conversion-factor dW --> dp/p
            z  = phase2z(conv_phase2z,pOUT,pOUTs)
            zp = W2zp(conv_W2zp,WOUT,dws)            
            
            DEBUG_MAP("(z,z')",(z,zp))

            T = T + dws       # soll-energy gained

            dbTab1Row = [      # for DEBUGGING
                '{:8.4f}'.format(WINs),
                '{:8.4g}'.format(WOUTs),
                '{:8.4f}'.format(r),
                '{:8.4f}'.format(Kr),
                '{:8.4f}'.format(i0),
                '{:8.4f}'.format(i1),
                '{:8.4g}'.format(WIN*1.e3),
                '{:8.4g}'.format(WOUT*1.e3),
                ]
            dbTab1Rows.append(dbTab1Row)
            dbTab2Row = [      # for DEBUGGING
                '{:8.4g}'.format(z*1.e3),
                '{:8.4g}'.format(zp*1.e2),
                '{:8.4g}'.format(xpf*1.e3),
                '{:8.4g}'.format(ypf*1.e3),
                '{:8.4g}'.format(degrees(pOUT-pOUTs)*1.e3),
                ]
            dbTab2Rows.append(dbTab2Row)

            # for DEBUGGING
            dbTab1Headr = ['WINs','WOUTs','r','Kr','i0','i1','WIN[KeV]','WOUT[KeV]']
            dbTab2Headr = ['z[mm]',"z'=dp/p[%]","x'[mrad]","y'[mrad]",'pOUT-pOUTs[mdeg]']
            DEBUG_MAP('\n'+(tblprnt(dbTab1Headr,dbTab1Rows)))
            DEBUG_MAP ('\n'+(tblprnt(dbTab2Headr,dbTab2Rows)))

            f_track = np.array([x,xpf,y,ypf,z,zp,T,1.,s,1.])
            return f_track
        # body --------------------------------------------------------------        
        for cnt,slice in enumerate(self.slices):
            DEBUG_MAP('SLICE # {} '.format(cnt),end='')  # for DEBUGGING
            f_track = slice_map(i_track)
            i_track = f_track
                # z = f_track[ZKOO]                # relativistic scaling. Is it needed?
                # betai = self.particle.beta
                # tkin  = self.particle.tkin
                # betaf = Proton(tkin=tkin+self.deltaW).beta
                # z = betaf/betai*z
                # f_track[ZKOO] = z
        return f_track
    def soll_map(self,i_track):
        def slice_map(track):
            """
            Mapping of soll track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
            """
            def Dwout_win(conv,i0,tk,sk,phi):
                # Formel 4.3.1 A.Shishlo
                return conv*i0*(tk*cos(phi)-sk*sin(phi))

            x        = track[XKOO]       # [0]
            xp       = track[XPKOO]      # [1]
            y        = track[YKOO]       # [2]
            yp       = track[YPKOO]      # [3]
            z        = track[ZKOO]       # [4] z
            zp       = track[ZPKOO]      # [5] dp/p
            T        = track[EKOO]       # [6] summe aller delta-energie
            s        = track[SKOO]       # [8] summe aller laengen

            particlei = slice.particle
            qV0       = slice.getV0()
            Tk        = slice.getTk()
            pINs      = slice.getPhase()   # soll phase @ (i)
            WINs      = slice.getW()       # soll energy @ (i)
            i0 = 1.

            dws = Dwout_win(qV0,i0,Tk,0.,pINs)    
            T = T + dws       # soll energy gained
            slice.setDW(dws)  # my delatW
            
            if (x,xp,y,yp,z,zp) != (0.,0.,0.,0.,0.,0.):
                DEBUG_SOLL_MAP('(T,dws,x,xp,y,yp,z,zp)',(T,dws,x,xp,y,yp,z,zp))
                raise RuntimeError("soll-track should have all zeros!")
                sys.exit(1)
            else:
                DEBUG_SOLL_MAP('(T,dws)',(T,dws))
            f_track = np.array([x,xp,y,yp,z,zp,T,1.,s,1.])
            return f_track
        # body --------------------------------------------------------------        
        for cnt,slice in enumerate(self.slices):
            DEBUG_SOLL_MAP('soll:SLICE # {} '.format(cnt),end='')  # for DEBUGGING
            f_track = slice_map(i_track)
            i_track = f_track
        return f_track
class TTFC(ELM._thin):
    # TODO 
    pass
def test0():
    from tracks import Track
    
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    
    ttfg = TTFG(length=0.044,gap=0.048,Ez=SF_tab)
    tkin = 50.
    ttfg.adapt_for_energy(tkin=tkin)
    # DEBUG_MODULE('ttfg.__dict__',ttfg.__dict__)      # for DEBUGGING
    slices = ttfg.slices
    for slice in slices:
        # DEBUG_MODULE('slice\n',slice.__dict__)      # for DEBUGGING
        pass

    z = 1.e-3
    x=y=1.e-2
    T = tkin
    # track-point fields:        x   x'  y  y'  z   z'  T  1   s   1
    track = Track(start=np.array([ x,  0., y, 0., z,  0., T, 1., 0., 1.]))
    ti = track.last()
    for i in range(1):
        # DEBUG_MAP(track.last_str())
        tf = ttfg.map(ti)
        track.append(tf)
        # DEBUG_MAP(track.last_str())
        ttfg.adapt_for_energy(tf[EKOO])    #enery adaptation
        ti = tf
if __name__ == '__main__':
    test0()