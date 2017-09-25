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
DEBUG_MODULE   = DEBUG_ON
DEBUG_MAP      = DEBUG_ON
DEBUG_SOLL_MAP = DEBUG_OFF
DEBUG_SLICE    = DEBUG_OFF
twopi          = 2*PI

## Transition Time Factors RF Gap Model
class TTFGslice(object):
    def __init__(self,ttfg,polyval,particle):
        self.elemnt     = ttfg    # the element this slice is part off
        self.dWf        = self.elemnt.dWf
        self.polyval    = polyval # polynom interval: ACHTUNG: E(z)=E0(1.+a*z+b*z**2), z in [cm] E0 in [MV/m]
        self.particle   = copy(particle)   # incoming sync. particle
        self.beta       = self.particle.beta
        self.gamma      = self.particle.gamma
        self.gb         = self.particle.gamma_beta
        self.k          = twopi/(self.elemnt.lamb*self.beta)
        self.Tk         = self._T (self.polyval,self.k)
        self.Tkp        = self._Tp(self.polyval,self.k)
        self.V0         = self._V0(self.polyval)
    def setSollPhase(self,value):
        self.phis = value
    def _T(self,poly,k):
        a  = poly.a
        b  = poly.b
        dz = poly.dz
        k  = k*1.e-2      # [1/m] --> [1/cm]
        f1 = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        t  = f1*f2
        # DEBUG_MAP('(T,k)',(t,k))
        return t
    def _Tp(self,poly,k):
        a   = poly.a
        b   = poly.b
        dz  = poly.dz
        k   = k*1.e-2      # [1/m] --> [1/cm]
        tp  = 2*sin(k*dz)/(k*(2*dz+2./3.*b*dz**3))
        tp  = tp*((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tp  = tp*1.e-2     # [cm] --> [m]
        return tp
    def _V0(self,poly):
        E0 = poly.E0                          #[MV/m]
        b  = poly.b                           #[1/cm**2]
        dz = poly.dz                          #[cm]
        v0 = (2*dz+2./3.*b*dz**3)*1.e-2       #[cm] --> [m]
        v0 = v0*E0*self.dWf
        return v0
    def mapSoll(self,i_track):
        i_track[EKOO] += self.deltaW
        return i_track
    def adjust_energy(self,tkin):
        self.particle(tkin)
        self.beta    = self.particle.beta
        self.gamma   = self.particle.gamma
        self.gb      = self.particle.gamma_beta
        self.k       = twopi/(self.elemnt.lamb*self.beta)
        self.Tk      = self._T (self.polyval,self.k)
        self.Tkp     = self._Tp(self.polyval,self.k)

        c     = PARAMS['lichtgeschwindigkeit']
        m0c3  = self.particle.e0*c
        omeg  = twopi*self.elemnt.freq
        I0    = 1.
        WIN   = self.particle.tkin
        DW    = self.V0*I0*self.Tk*cos(self.phis)
        WOUT  = WIN+DW
        PHIN  = self.phis
        DPHI  = self.V0*omeg/(m0c3*self.gb**3)*I0*self.Tk*sin(self.phis)
        PHOUT = PHIN+DPHI

        self.WIN      = WIN
        self.WOUT     = WOUT
        self.deltaW   = DW
        self.PHIN     = PHIN
        self.PHOUT    = PHOUT
        self.deltaPHI = DPHI
        return 
    def wout_minus_win(self,conv,i0,tk,sk,phi):
        # Formel 4.3.1 A.Shishlo
        return conv*i0*(tk*cos(phi)-sk*sin(phi))
    def phiout_minus_phiin(self,fac,gamma,r,i0,i1,tk,sk,tkp,skp,phi):
        # Formel 4.3.2 A.Shishlo
        return  fac*i0*(tkp*sin(phi)+skp*cos(phi)+gamma*r*i1*(tk*sin(phi)+sk*cos(phi)))
    def map(self,i_track):
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z~(phi-phis)
        zp       = i_track[ZPKOO]      # [5] dp/p~dT
        T        = i_track[EKOO]       # [6] summe aller dT
        s        = i_track[SKOO]       # [8] summe aller laengen
        c          = PARAMS['lichtgeschwindigkeit']
        m0c2       = self.particle.e0
        m0c3       = m0c2*c
        omeg       = twopi*self.elemnt.freq
        be         = self.particle.beta
        ga         = self.particle.gamma
        gb         = self.particle.gamma_beta
        k          = omeg/(c*be)
        Tk         = self._T (self.polyval,k)
        Tkp        = self._Tp(self.polyval,k)
        r          = sqrt(x**2+y**2)                # radial coordinate
        Kbess      = omeg/(c*gb) * r
        i0         = I0(Kbess)                      # bessel function
        i1         = I1(Kbess)                      # bessel function
        fact       = self.V0*omeg/(m0c3*gb**3)

        pin       = -z*omeg/(c*be) + self.PHIN     # phase @ (i)  ~ (-z)
        win       = (zp*(ga+1.)/ga+1.)*self.WIN    # energy @ (i) ~ (zp)
        
        wowi = self.wout_minus_win(self.V0,i0,Tk,0.,pin) 
        wout = win + wowi                              # energy @ (f)

        popi = self.phiout_minus_phiin(fact,ga,r,i0,i1,Tk,0.,Tkp,0.,pin)
        pout  = pin + popi                             # phase @ (f)
        
        dp = +(pout-self.PHOUT)    # delta phase @ (f)
        dw = +(wout-self.WOUT)    # delta energy @ (f)
        
        zf  = -dp*(c*be)/omeg
        zpf = ga/(ga+1)*dw/self.WOUT

        f_particle = copy(self.particle)(tkin=wout)
        gbf        = f_particle.gamma_beta

        fact = self.V0/(m0c2*gb*gbf)*i1
        gbi = self.gb
        if r > 0.:
            xp = gbi/gbf*xp-x/r*fact*Tk*sin(pin)
            yp = gbi/gbf*yp-y/r*fact*Tk*sin(pin)
        elif r == 0.:
            xp = gbi/gbf*xp
            yp = gbi/gbf*yp

        T = T + self.deltaW       # soll-energy gained
        f_track = np.array([x,xp,y,yp,zf,zpf,T,1.,s,1.])

        dbTab1Row = [      # for DEBUGGING
            '{:8.4f}'.format(degrees(pout)),
            '{:8.4f}'.format(degrees(pin)),
            '{:8.4g}'.format(degrees(popi)),
            '{:8.4g}'.format(dp),
            '{:8.4g}'.format(wout),
            '{:8.4g}'.format(win),
            '{:8.4g}'.format(wowi),
            '{:8.4g}'.format(self.WOUT),
            '{:8.4g}'.format(dw),
            '{:8.3f}'.format(self.V0*1.e3),
            ]
        self.elemnt.dbTab1Rows.append(dbTab1Row)

        #   dbTab2Row = [      # for DEBUGGING
        #     '{:8.4g}'.format(x*1.e3),
        #     '{:8.4g}'.format(xp*1.e2),
        #     '{:8.4g}'.format(y*1.e3),
        #     '{:8.4g}'.format(yp*1.e3),
        #     '{:8.4g}'.format(z*1.e3),
        #     '{:8.4g}'.format(zp*1.e3),
        #     '{:8.4g}'.format(r*1.e3),
        #     '{:8.4g}'.format(TINs),
        #     '{:8.3g}'.format(TpINs),
        #     '{:8.4g}'.format(i0-1.),
        #     '{:8.4g}'.format(i1),
        #     ]
        # dbTab2Rows.append(dbTab2Row)
        return f_track
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
            raise RuntimeError("can't create object without SF-data! - STOP")
            sys.exit(1)
    def adjust_energy(self,tkin):
        def adjust_slice_energy(tkin):
            next_phase = self.phis
            next_tkin  = tkin
            Tklist = []
            for slice in self.slices:
                slice.setSollPhase(next_phase)
                slice.adjust_energy(next_tkin)
                Tklist.append(slice.Tk)
                next_phase = slice.PHOUT
                next_tkin  = slice.WOUT
            deltaW   = next_tkin-tkin
            self.tr  = min(Tklist)
            return deltaW
        # body --------------------------------------------------------------        
        self.particle(tkin)                     # adjust my particle energy
        DEBUG_SLICE('adjust_energy       PARTICLE      :',self.particle)
        DEBUG_SLICE('adjust_energy       PARTICLE tkin :',self.particle.tkin)
        self.deltaW = adjust_slice_energy(tkin) # ajust slices
        self.matrix[EKOO,DEKOO] = self.deltaW   # my deltaW to lin.map
        return(self)
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
    def make_slices(self,anz=0):
        res = [self]
        DEBUG_SLICE('ttfg-make_slices: ',res)
        return res
    def shorten(self,l=0.):
        return self
    def map(self,i_track):
        """
        Mapping of track from position (i) to (f)
        """
        if FLAGS['map']:
            # NOTE: mapping with RFB-map
            f_track = self._map(i_track)
        else:
            # NOTE: linear mapping with T3D matrix
            f_track = super().map(i_track)
        return f_track
    def _map(self,i_track):
        # def slice_map(track,index):
        #     """
        #     Mapping of track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
        #     """
        #     # def phiout_minus_phiin(conv,g,r,i0,i1,tkp,skp,tk,sk,phi):
        #     #     # Formel 4.3.2 A.Shishlo
        #     #     return  conv*i0*(tkp*sin(phi)+skp*cos(phi)+g*r*i1*(tk*sin(phi)+sk*cos(phi)))
        #     # def wout_minus_win(conv,i0,tk,sk,phi):
        #     #     # Formel 4.3.1 A.Shishlo
        #     #     return conv*i0*(tk*cos(phi)-sk*sin(phi))

       #       x        = track[XKOO]       # [0]
        #     xp       = track[XPKOO]      # [1]
        #     y        = track[YKOO]       # [2]
        #     yp       = track[YPKOO]      # [3]
        #     z        = track[ZKOO]       # [4] z~(phi-phis)
        #     zp       = track[ZPKOO]      # [5] dp/p~dT
        #     T        = track[EKOO]       # [6] summe aller dT
        #     s        = track[SKOO]       # [8] summe aller laengen
        #     
        #     # DEBUG_MAP('(z,zp)',(z,zp))

       #       c          = PARAMS['lichtgeschwindigkeit']
        #     m0c2       = self.particle.e0
        #     m0c3       = m0c2*c
        #     omeg       = twopi*self.freq
        #     
        #     soll = slice.particle                  # soll @ (i)
        #     besi       = soll.beta
        #     gasi       = soll.gamma
        #     gbsi       = soll.gamma_beta
        #     qV0       = slice.getV0()
        #     ksi        = omeg/(c*besi)
        #     TINs      = slice.getT(ksi)
        #     TpINs     = slice.getTp(ksi)
        #     PIN       = slice.getPhaseIN()           # soll phase  @ (i)
        #     WIN       = slice.getWIN()               # soll energy @ (i)

       #       pin       = -z*omeg/(c*besi) + PIN     # phase @ (i)  ~ (-z)
        #     win       = (zp*(gasi+1.)/gasi+1.)*WIN # energy @ (i) ~ (zp)
        #     
        #     # DEBUG_MAP('win ',win)
        #     
        #     particle  = copy(soll)(tkin=win)       # particle @ (i)
        #     bei        = particle.beta
        #     gai        = particle.gamma
        #     gbi        = particle.gamma_beta
        #     ki         = omeg/(c*bei)
        #     TIN       = slice.getT(ki)
        #     TpIN      = slice.getTp(ki)
        #     
        #     r = sqrt(x**2+y**2)                 # radial coordinate
        #     Kbess = omeg/(c*gbi) * r
        #     i0 = I0(Kbess)                      # bessel function
        #     i1 = I1(Kbess)                      # bessel function

       #       WOUT = WIN + wout_minus_win(qV0,1.,TINs,0.,PIN)    # soll energy @ (f)
        #     wout = win + wout_minus_win(qV0,i0,TIN ,0.,pin)    # energy @ (f)

       #       # DEBUG_MAP('(WIN,WOUT,WIN+DW,wout)',(WIN,WOUT,WIN+slice.getDW(),wout))

       #       facs  = qV0*omeg/(m0c3*gbsi**3)
        #     fac   = qV0*omeg/(m0c3*gbi*3)
        #     POUT  = PIN + phiout_minus_phiin(facs,gasi ,0.,1.,0.,TINs,0.,TpINs,0.,PIN)  # soll phase  @ (f)
        #     if index == len(self.slices)-1:
        #         pass
        #     else:
        #         self.slices[index+1].setPhase(POUT)
        #     pout  = pin + phiout_minus_phiin(facs ,gai  ,r ,i0,i1,TIN ,0.,TpIN ,0.,pin)  # phase @ (f)
        #     
        #     DEBUG_MAP('(PIN,POUT,pout) ',(PIN,POUT,pout))

       #       dp = +(pout-POUT)    # delta phase @ (f)
        #     dw = +(wout-WOUT)    # delta energy @ (f)
        #     
        #     zf  = -dp*(c*bei)/omeg
        #     zpf = gai/(gai+1)*dw/WOUT

       #       particle  = copy(soll)(tkin=wout)
        #     gbf       = particle.gamma_beta

       #       common = qV0/(m0c2*gbi*gbf)*i1
        #     if r > 0.:
        #         xp = gbi/gbf*xp-x/r*common*TIN*sin(pin)
        #         yp = gbi/gbf*yp-y/r*common*TIN*sin(pin)
        #     elif r == 0.:
        #         xp = gbi/gbf*xp
        #         yp = gbi/gbf*yp

       #       T = T + WOUT-WIN       # soll-energy gained
        #     f_track = np.array([x,xp,y,yp,zf,zpf,T,1.,s,1.])
        #     
        #     # DEBUG utilities
        #     qV0_av.append(qV0) # gather gap voltages for average

       #       dbTab1Row = [      # for DEBUGGING
        #         '{:8.4f}'.format(degrees(pout)),
        #         '{:8.4f}'.format(degrees(pin)),
        #         '{:8.4g}'.format(degrees(pout-pin)),
        #         '{:8.4g}'.format(dp),
        #         '{:8.4g}'.format(wout),
        #         '{:8.4g}'.format(win),
        #         '{:8.4g}'.format(wout-win),
        #         '{:8.4g}'.format(WOUT),
        #         '{:8.4g}'.format(dw),
        #         '{:8.3f}'.format(qV0*1.e3),
        #         ]
        #     dbTab1Rows.append(dbTab1Row)

       #       dbTab2Row = [      # for DEBUGGING
        #         '{:8.4g}'.format(x*1.e3),
        #         '{:8.4g}'.format(xp*1.e2),
        #         '{:8.4g}'.format(y*1.e3),
        #         '{:8.4g}'.format(yp*1.e3),
        #         '{:8.4g}'.format(z*1.e3),
        #         '{:8.4g}'.format(zp*1.e3),
        #         '{:8.4g}'.format(r*1.e3),
        #         '{:8.4g}'.format(TINs),
        #         '{:8.3g}'.format(TpINs),
        #         '{:8.4g}'.format(i0-1.),
        #         '{:8.4g}'.format(i1),
        #         ]
        #     dbTab2Rows.append(dbTab2Row)
        #     return f_track
        # body --------------------------------------------------------------        
        # table #1
        self.dbTab1Rows  = []      # for DEBUGGING
        self.dbTab1Headr = ['pout','pin','pout-pin','dp=pout-POUT','wout','win','wout-win','WOUT','dw=wout-WOUT','qV0*10^3']
        # table #2
        # dbTab2Rows  = []      # for DEBUGGING
        # dbTab2Headr = ['x*10^3','xp*10^3','y*10^3','yp*10^3','z*10^3','zp*10^3','r*10^3','TINs',"TpINs",'i0-1','i1']      # for DEBUGGING

        for cnt,slice in enumerate(self.slices):
            # DEBUG_MAP('ttfg-map: {} tkin {} '.format(cnt,self.particle.tkin),slice)
            f_track = slice.map(i_track)
            i_track = f_track
                # z = f_track[ZKOO]                # relativistic scaling. Is it needed?
                # betai = self.particle.beta
                # tkin  = self.particle.tkin
                # betaf = Proton(tkin=tkin+self.deltaW).beta
                # z = betaf/betai*z
                # f_track[ZKOO] = z
        DEBUG_MAP('ttfg-map: settings of slices')
        # DEBUG_MAP('<qV0>/slice*10^3 {:8.4f}\n'.format(qV0_av),tblprnt(dbTab1Headr,dbTab1Rows))
        DEBUG_MAP('\n',tblprnt(self.dbTab1Headr,self.dbTab1Rows))
        # DEBUG_MAP('\n',tblprnt(dbTab2Headr,dbTab2Rows))
        return f_track
    def soll_map(self,i_track):
        # def slice_map(track):
        #     """
        #     Mapping of soll track from position (i) to (f) in TTF-Gap model approx. (A.Shislo 4.4)
        #     """
        #     def wout_minus_win(conv,i0,tk,sk,phi):
        #         # Formel 4.3.1 A.Shishlo
        #         return conv*i0*(tk*cos(phi)-sk*sin(phi))

       #       x        = track[XKOO]       # [0]
        #     xp       = track[XPKOO]      # [1]
        #     y        = track[YKOO]       # [2]
        #     yp       = track[YPKOO]      # [3]
        #     z        = track[ZKOO]       # [4] z
        #     zp       = track[ZPKOO]      # [5] dp/p
        #     T        = track[EKOO]       # [6] summe aller delta-energie
        #     s        = track[SKOO]       # [8] summe aller laengen

       #       qV0       = slice.getV0()
        #     Tk        = slice.getTk()
        #     PIN       = slice.getPhaseIN()   # soll phase @ (i)

       #       dws = wout_minus_win(qV0,1.,Tk,0.,PIN)    
        #     T = T + dws       # soll energy gained
        #     slice.setDW(dws)  # my deltaW
        #     
        #     if (x,xp,y,yp,z,zp) != (0.,0.,0.,0.,0.,0.):
        #         DEBUG_SOLL_MAP('(T,dws,x,xp,y,yp,z,zp)',(T,dws,x,xp,y,yp,z,zp))
        #         raise RuntimeError("soll-track should have all zeros! - STOP")
        #         sys.exit(1)
        #     else:
        #         DEBUG_SOLL_MAP('(T,dws)',(T,dws))
        #     f_track = np.array([x,xp,y,yp,z,zp,T,1.,s,1.])
        #     return f_track
        # body --------------------------------------------------------------        
        for cnt,slice in enumerate(self.slices):
            DEBUG_SOLL_MAP('soll:SLICE # {} '.format(cnt),end='')  # for DEBUGGING
            f_track = slice.mapSoll(i_track)
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
    
    ttfg = TTFG(gap=0.048,Ez=SF_tab)
    tkin = 50.
    ttfg.adjust_energy(tkin=tkin)
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
        ttfg.adjust_energy(tf[EKOO])    #enery adaptation
        ti = tf
def test1():
    from tracks import Track
    
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    
    ttfg = TTFG(PhiSoll=radians(-20.),gap=0.048,Ez=SF_tab)
    tkin = 150.
    ttfg.adjust_energy(tkin=tkin)
    DEBUG_MODULE('ttfg.__dict__',ttfg.__dict__)      # for DEBUGGING
    slices = ttfg.slices
    
    von = 0.
    bis = 1.
    anz = 20
    delta = (bis-von)/anz
    x=xp=y=yp=z=zp=0.0
    x = 1.e-1
    z = von
    T = 0.
    # start            x   x'  y  y'  z   z'  T  1   s   1
    start = np.array([ x,  xp, y, yp, z,  zp, T, 1., 0., 1.])
    trck = Track(start=start)
    ti = trck.last()
    for i in range(anz+1):
        # DEBUG_MODULE(trck.last_str())
        tf = ttfg.map(ti)
        trck.append(tf)
        # print(trck.last_str())
        z += delta
        start[4] = z
        ti = start
    print(trck.asTable())
if __name__ == '__main__':
    # test0()    
    test1()