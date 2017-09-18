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
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
def DEBUG_MAP_ON(*args):
    DEBUG(*args)
def DEBUG_MAP_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF
DEBUG_MAP    = DEBUG_MAP_OFF
twopi        = 2*PI

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
s    def _T(self,poly,k):
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
            for slice in self.slices:
                slice.setPhase(next_phase)
                slice.adapt_for_energy(next_tkin)
                DPhase     = slice.getDPhase()
                DW         = slice.getDW()
                next_phase = slice.getPhase()+DPhase
                next_tkin  = slice.particle.tkin+DW
        
        self.__init__(PhiSoll     = self.phis,
                      fRF         = self.freq,
                      label       = self.label,
                      particle    = self.particle(tkin),
                      gap         = self.gap,
                      Ez          = self.Ez,
                      dWf         = self.dWf)

        adapt_slices_for_energy(tkin)
        deltaW = self.slices[-1].getW()-self.slices[0].getW()  #deltaW for lin.map
        self.matrix[EKOO,DEKOO] = deltaW
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
            # zur Erinnerung:
            #     sync.Teilchen: (W0,DW0,ph0) = (tkin,delta(tkin),sync.phase)
            #     Teilchen: ( W  ,  ph  , z  <-->  Dph ,  zp <-->  DW ,       DDph    ,  DDW      )=
            #               (Tkin,phase,z-z0<-->phi-ph0, dp/p<--> W-W0, delta(phi-ph0),delta(W-W0))
            x        = track[XKOO]       # [0]
            xp       = track[XPKOO]      # [1]
            y        = track[YKOO]       # [2]
            yp       = track[YPKOO]      # [3]
            z        = track[ZKOO]       # [4] z
            zp       = track[ZPKOO]      # [5] dp/p
            T        = track[EKOO]       # [6] summe aller delta-energie
            s        = track[SKOO]       # [8] summe aller laengen

            DEBUG_MAP("z,z' ",(z,zp))

            m0c2       = self.particle.e0
            lamb       = self.lamb
            c          = PARAMS['lichtgeschwindigkeit']

            DEBUG_MAP('(m0c2,lamb,c)',(m0c2,lamb,c))

            dbTab1Rows  = []      # for DEBUGGING
            dbTab2Rows  = []      # for DEBUGGING

            particlei = slice.particle
            betai     = particlei.beta
            gammai    = particlei.gamma
            gbi       = particlei.gamma_beta
            qV0       = slice.getV0()
            Tk        = slice.getTk()
            Tkp       = slice.getTkp()
            Sk        = 0.
            Skp       = 0.

            DEBUG_MAP('(betai,gammai,gbi,qV0,Tk,Tkp)',(betai,gammai,gbi,qV0,Tk,Tkp))

            r = sqrt(x**2+y**2)              # radial coordinate
            Kr = (twopi*r)/(lamb*gbi)
            i0 = I0(Kr)                      # bessel function
            i1 = I1(Kr)                      # bessel function

            DEBUG_MAP('(r,Kr,i0,i1)',(r,Kr,i0,i1))

            conv_zDPhi  = -twopi/(betai*lamb)       # conversion factor z --> Dph @ (in)
            conv_dPdT   = (m0c2*betai**2*gammai)    # conversion factor dp/p --> dT @ (in)

            DEBUG_MAP('(conv_zDPhi,conv_dPdT)',(conv_zDPhi,conv_dPdT))

            Dph   = conv_zDPhi*z            # phase(in): z --> dphi
            ph0   = slice.getPhase()        # phase(in): ph0

            DW    = conv_dPdT*zp            # W(in): dp/p --> dT
            W0    = slice.getW()            # W(in): W0

            DEBUG_MAP('(Dph,DW,ph0,W0)',(Dph,DW,ph0,W0))

            commonf = twopi*qV0/(lamb*m0c2*gbi**3)
            ph = Dph+ph0
            DDph=phfmphi = commonf*(i0*(Tkp*sin(ph)+Skp*cos(ph))+
                gammai*r*i1*(Tk*sin(ph)+Sk*cos(ph)))            # Formel 4.3.2 A.Shishlo

            DEBUG_MAP('(commonf,ph,phfmphi)',(commonf,ph,phfmphi))

            W = W0 + DW                     #particle tkin
            if W < 0.:
                raise RuntimeError('negative kinetic energy {:8.4g}'.format(Wi))
                sys.exit(1)
            DW0 = slice.getDW()
            W0 = W0+DW0                     # soll tkin
            particlef = copy(particlei)(tkin=W0)
            betaf = particlef.beta
            gammaf = particlef.gamma
            gbf = particlef.gamma_beta

            DDW=WfmWi= qV0*i0*(Tk*cos(ph)-Sk*sin(ph))  # Formel 4.3.1 A.Shishlo

            DEBUG_MAP('(W,DW0,W0,betaf,gammaf,gbf,WfmWi)',(W,DW0,W0,betaf,gammaf,gbf,WfmWi))

            conv_Dphiz = -(betaf*lamb)/twopi        # conversion factor Dph --> z @ (f)
            conv_dTdP  = 1./(m0c2*betaf**2*gammaf)  # conversion factor dT --> dp/p @ (f)
            z  = conv_Dphiz*(DDph+Dph)              # phase(out): dT --> dp/p
            zp = conv_dTdP* (DDW+DW)                # W(out): Dph --> z

            # z = gbf/gbi*z    relativistic scaling need?
            
            DEBUG_MAP("(conv_Dphiz,conv_dTdP)",(conv_Dphiz,conv_dTdP))
            DEBUG_MAP("(z,z')",(z,zp))

            commonf = qV0/(m0c2*gbi*gbf)*i1         # common factor
            if r > 0.:
                xpf = gbi/gbf*xp-x/r*commonf*(Tk*sin(ph)+Sk*cos(ph)) # Formel 4.3.3 A.Shishlo
                ypf = gbi/gbf*yp-y/r*commonf*(Tk*sin(ph)+Sk*cos(ph))
            elif r == 0.:
                xpf = gbi/gbf*xp
                ypf = gbi/gbf*yp

            T = T + DW0       # soll energy gained

            dbTab1Row = [      # for DEBUGGING
                '{:8.4f}'.format(W0),
                '{:8.4g}'.format((DW0)*1.e3),
                '{:8.4f}'.format(r),
                '{:8.4f}'.format(Kr),
                '{:8.4f}'.format(i0),
                '{:8.4f}'.format(i1),
                '{:8.4g}'.format(DW*1.e3),
                '{:8.4g}'.format(degrees(Dph)),
                ]
            dbTab1Rows.append(dbTab1Row)
            dbTab2Row = [      # for DEBUGGING
                '{:8.4g}'.format(z*1.e3),
                '{:8.4g}'.format(zp*1.e2),
                '{:8.4g}'.format(xpf*1.e3),
                '{:8.4g}'.format(ypf*1.e3),
                '{:8.4g}'.format(WfmWi*1.e3),
                '{:8.4g}'.format(degrees(phfmphi)*1.e3),
                ]
            dbTab2Rows.append(dbTab2Row)

            # for DEBUGGING
            dbTab1Headr = ['W0','DW0[KeV]','r','Kr','i0','i1','DW[KeV]','Dph[deg]']
            dbTab2Headr = ['z[mm]',"z'=dp/p[%]","x'[mrad]","y'[mrad]",'DDW[KeV]','DDphi[mdeg]']
            DEBUG_MODULE('at position IN:\n'+(tblprnt(dbTab1Headr,dbTab1Rows)))
            DEBUG_MODULE('at position OUT:\n'+(tblprnt(dbTab2Headr,dbTab2Rows)))

            f_track = np.array([x,xpf,y,ypf,z,zp,T,1.,s,1.])
            return f_track

        for cnt,slice in enumerate(self.slices):
            DEBUG_MAP('\nSLICE # {}'.format(cnt))  # for DEBUGGING
            f_track = slice_map(i_track)
            i_track = f_track
        return f_track
# class TTFC(ELM._thin):
#     """
#     Transition Time Factors RF Gap Model (A.Shishlo ORNL/TM-2015/247)
#     MAP: D.map()*TTFG.map()*D.map()
#     """
#     def __init__(self,
#                     PhiSoll    = radians(PARAMS['soll_phase']),
#                     fRF        = PARAMS['frequenz'],
#                     label      = 'TTFG',
#                     particle   = PARAMS['sollteilchen'],
#                     gap        = PARAMS['spalt_laenge'],
#                     length     = 0.,
#                     Ez         = None,                # return of SFdata
#                     dWf        = FLAGS['dWf']):
#         super().__init__(particle=particle)
#         if Ez == None: raise RuntimeError('TTFG: missing E(z) table')
#         self.viseo  = 0.25
#         self.phis   = PhiSoll    # [radians] soll phase
#         self.freq   = fRF        # [Hz]  RF frequenz
#         self.label  = label
#         self.gap    = gap        # [m]
#         self.length = length     # [m]
#         self.dWf    = dWf
#         self.lamb   = PARAMS['wellenlänge']
#         self.Ez     = Ez         # the SF-table, a list(DataPoints)
#         try:
#             self.Epeak  = Ez.Epeak
#             self.slices = self._make_slices()
#         except:
#             raise RuntimeError("can't create object without SF-data!")
#             sys.exit(1)
#         
#         # drifts before and after gap
#         self.inD     = ELM.D(length=0.5*length,label='dtI', particle=particle) 
#         self.outD    = ELM.D(length=0.5*length,label='dtF', particle=particle)
# 
#         self.matrix  = ELM.D(length=length,label='ttfg',particle=particle).() 
#     def adapt_for_enery(self,tkin):
#         def adapt_slices_for_energy(self,tkin):
#             next_phase = self.phis
#             next_tkin  = tkin
#             for cnt,slice in enumerate(self.slices):
#                 slice.setPhase(next_phase)
#                 slice.adapt_for_energy(next_tkin)
#                 DPhase     = slice.getDPhase()
#                 DW         = slice.getDW()
#                 next_phase = slice.getPhase()+DPhase
#                 next_tkin  = slice.particle.tkin+DW
#         
#         self.__init__(PhiSoll     = self.phis,
#                       fRF         = self.freq,
#                       label       = self.label,
#                       particle    = self.particle(tkin),
#                       gap         = self.gap,length=self.
#                       length      = self.length,
#                       Ez          = self.Ez,
#                       dWf         = self.dWf)
#         self.adapt_slices_for_energy(tkin)
#         # self.deltaW is delta-W of sync.particle through TTF-gap
#         self.deltaW = self.slices[-1].particle.tkin-self.slices[0].particle.tkin
#         return(self)
#     def _make_slices(self,anz=0):
#         slices = []
#         zl = -self.gap/2.*100.   # [m] --> [cm]
#         zr = -zl
#         for poly in self.Ez.Ez_poly():
#             zil = poly.zl
#             zir = poly.zr
#             if zil < zl or zir > zr: continue
#             slice = TTFGslice(self,poly,self.particle)  # instanciate TTFGslices
#             slices.append(slice)
#         return slices
#     def make_slices(self,anz=0):
#         return [self.inD,self,self.outD]
#     def map(self,i_track):
#         def slice_map(track):
#             """
#             Mapping of track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
#             """
#             x        = track[XKOO]       # [0]
#             xp       = track[XPKOO]      # [1]
#             y        = track[YKOO]       # [2]
#             yp       = track[YPKOO]      # [3]
#             z        = track[ZKOO]       # [4] z
#             zp       = track[ZPKOO]      # [5] dp/p
#             T        = track[EKOO]       # [6] summe aller deltaW
#             s        = track[SKOO]       # [8] summe aller laengen
# 
#             DEBUG_MAP("z,z' ",(z,zp))
# 
#             m0c2       = self.particle.e0
#             lamb       = self.lamb
#             c          = PARAMS['lichtgeschwindigkeit']
# 
#             DEBUG_MAP('(m0c2,lamb,c)',(m0c2,lamb,c))
# 
#             dbTab1Rows  = []      # for DEBUGGING
#             dbTab2Rows  = []      # for DEBUGGING
# 
#             particlei = slice.particle
#             betai     = particlei.beta
#             gammai    = particlei.gamma
#             gbi       = particlei.gamma_beta
#             qV0       = slice.getV0()
#             Tk       = slice.getTk()
#             Tkp      = slice.getTkp()
#             Sk       = 0.
#             Skp      = 0.
# 
#             DEBUG_MAP('(betai,gammai,gbi,qV0,Tk,Tkp)',(betai,gammai,gbi,qV0,Tk,Tkp))
# 
#             r = sqrt(x**2+y**2)              # radial coordinate
#             Kr = (twopi*r)/(lamb*gbi)
#             i0 = I0(Kr)                      # bessel function
#             i1 = I1(Kr)                      # bessel function
# 
#             DEBUG_MAP('(r,Kr,i0,i1)',(r,Kr,i0,i1))
# 
#             conv_zDPhi  = -twopi/(betai*lamb)       # conversion factor z --> Dph @ (in)
#             conv_dPdT   = (m0c2*betai**2*gammai)    # conversion factor dp/p --> dT @ (in)
# 
#             DEBUG_MAP('(conv_zDPhi,conv_dPdT)',(conv_zDPhi,conv_dPdT))
# 
#             Dph   = conv_zDPhi*z            # phase(in): z --> dphi
#             ph0   = slice.getPhase()        # phase(in): ph0
# 
#             DW    = conv_dPdT*zp            # W(in): dp/p --> dT
#             W0    = slice.getW()            # W(in): W0
# 
#             DEBUG_MAP('(Dph,DW,ph0,W0)',(Dph,DW,ph0,W0))
# 
#             commonf = twopi*qV0/(lamb*m0c2*gbi**3)
#             ph = Dph+ph0
#             phfmphi = commonf*(i0*(Tkp*sin(ph)+Skp*cos(ph))+
#                 gammai*r*i1*(Tk*sin(ph)+Sk*cos(ph)))            # Formel 4.3.2 A.Shishlo
# 
#             DEBUG_MAP('(commonf,ph,phfmphi)',(commonf,ph,phfmphi))
# 
#             W = W0 + DW                     # W(in)
#             if W < 0.:
#                 raise RuntimeError('negative kinetic energy {:8.4g}'.format(Wi))
#                 sys.exit(1)
#             DW0 = slice.getDW()
#             W0f = W0+DW0                    # W(out): W0
#             particlef = copy(particlei)(tkin=W0f)
#             betaf = particlef.beta
#             gammaf = particlef.gamma
#             gbf = particlef.gamma_beta
# 
#             WfmWi= qV0*i0*(Tk*cos(ph)-Sk*sin(ph))  # Formel 4.3.1 A.Shishlo
# 
#             DEBUG_MAP('(W,DW0,W0f,betaf,gammaf,gbf,WfmWi)',(W,DW0,W0f,betaf,gammaf,gbf,WfmWi))
# 
#             conv_Dphiz = -(betaf*lamb)/twopi        # conversion factor Dph --> z @ (f)
#             conv_dTdP  = 1./(m0c2*betaf**2*gammaf)  # conversion factor dT --> dp/p @ (f)
#             z  = conv_Dphiz*(phfmphi+Dph)           # phase(out): dT --> dp/p
#             zp = conv_dTdP*(WfmWi+DW)               # W(out): Dph --> z
# 
#             DEBUG_MAP("(conv_Dphiz,conv_dTdP)",(conv_Dphiz,conv_dTdP))
#             DEBUG_MAP("(z,z')",(z,zp))
# 
#             commonf = qV0/(m0c2*gbi*gbf)*i1         # common factor
#             if r > 0.:
#                 xpf = gbi/gbf*xp-x/r*commonf*(Tk*sin(ph)+Sk*cos(ph)) # Formel 4.3.3 A.Shishlo
#                 ypf = gbi/gbf*yp-y/r*commonf*(Tk*sin(ph)+Sk*cos(ph))
#             elif r == 0.:
#                 xpf = gbi/gbf*xp
#                 ypf = gbi/gbf*yp
# 
#             T = T + DW0       # W0: energy gained
# 
#             dbTab1Row = [      # for DEBUGGING
#                 '{:8.4f}'.format(W0),
#                 '{:8.4g}'.format((DW0)*1.e3),
#                 '{:8.4f}'.format(r),
#                 '{:8.4f}'.format(Kr),
#                 '{:8.4f}'.format(i0),
#                 '{:8.4f}'.format(i1),
#                 '{:8.4g}'.format(DW*1.e3),
#                 '{:8.4g}'.format(degrees(Dph)),
#                 ]
#             dbTab1Rows.append(dbTab1Row)
#             dbTab2Row = [      # for DEBUGGING
#                 '{:8.4g}'.format(z*1.e3),
#                 '{:8.4g}'.format(zp*1.e2),
#                 '{:8.4g}'.format(xpf*1.e3),
#                 '{:8.4g}'.format(ypf*1.e3),
#                 '{:8.4g}'.format(WfmWi*1.e3),
#                 '{:8.4g}'.format(degrees(phfmphi)*1.e3),
#                 ]
#             dbTab2Rows.append(dbTab2Row)
# 
#             # for DEBUGGING
#             dbTab1Headr = ['W0','DW0[KeV]','r','Kr','i0','i1','DW[KeV]','Dph[deg]']
#             dbTab2Headr = ['z[mm]',"z'=dp/p[%]","x'[mrad]","y'[mrad]",'DDW[KeV]','DDphi[mdeg]']
#             DEBUG_MODULE('at position IN:\n'+(tblprnt(dbTab1Headr,dbTab1Rows)))
#             DEBUG_MODULE('at position OUT:\n'+(tblprnt(dbTab2Headr,dbTab2Rows)))
# 
#             f_track = np.array([x,xpf,y,ypf,z,zp,T,1.,s,1.])
#             return f_track
# 
#         for cnt,slice in enumerate(self.slices):
#             DEBUG_MAP('\nSLICE # {}'.format(cnt))  # for DEBUGGING
#             f_track = slice_map(i_track)
#             i_track = f_track
#         f_track[SKOO] += self.length         # ACHTUNG: add TTGF's length!
#         return f_track
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