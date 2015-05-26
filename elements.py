#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import wille,CONF,dictprnt,objprnt,Beam,k0,dBdz,scalek0,printv
import numpy as NP 
from math import sqrt,sinh,cosh,sin,cos,fabs,tan,floor,modf,pi,radians
from copy import copy

class _matrix(object): ## the mother of all 6x6 matrices
    _dim = 6   # 6x6 matrices
    def __init__(self):
        self.matrix=NP.eye(_matrix._dim)    ## 6x6 unit matrix
        self.label=''
        self.length=0.         ## default zero length!
        self.slice_min = 0.01  ## minimal slice length
        self.viseo = 0.
    def out(self):
        printv(1,self.label)
        printv(1,self.matrix)
    def __mul__(self,other):
        product=NP.einsum('ij,jk',self.matrix,other.matrix)
        res=_matrix()
        if (self.label == ''):
            res.label=other.label
        else:
            res.label=self.label+'*'+other.label
        res.length=self.length+other.length
        res.matrix=product
        return res
    def reverse(self):
        raise RuntimeError('_matrix:reverse not ready yet!')
        res=_matrix()
        for i in range(_matrix._dim):
            for k in range(_matrix._dim):
                res.matrix[i][k] = self.matrix[i][k]
        res.matrix[0][0] = self.matrix[1][1]
        res.matrix[1][1] = self.matrix[0][0]
        res.matrix[2][2] = self.matrix[3][3]
        res.matrix[3][3] = self.matrix[2][2]
        res.label = '('+self.label+')r'
        return res
    def trace(self):
        return self.tracex()+self.tracey
    def tracex(self):
        res = 0.
        for i in range(2):
            res += self.matrix[i,i]
        return res
    def tracey(self):
        res = 0.
        for i in range(2,4):
            res += self.matrix[i,i]
        return res
    def shorten(self,length=0.):
        return self
    def step_through(self,anz=10):
        """
        Step through an element - the central nontrivial function.
        Default is 10 steps/element.
        Minimal step size is self.slice_min.
        """
        step = self.length/anz
        if step < self.slice_min:
            step  = self.slice_min
        if self.length == 0.:           ## zero length element (like WD or CAV)
            anz = 0; fanz= 0.; rest= 0.
            mr=self
        else:
#            fanz = self.length/step
#            anz  = floor(fanz)
#            rest = self.length - anz*step
            (rest,fanz) = modf(self.length/step)
            anz = int(fanz)
            rest = self.length * rest
            mx = self.shorten(step)
            if fabs(rest) > 1.e-9:
                mr = self.shorten(rest)
            else:
                mr=I(label=self.label,viseo=self.viseo)
        # print('label={} fanz={} anz={} step={} rest={}'.format(self.label,fanz,anz,step,rest)) 
        for i in range(anz+1):
            if i == anz:
                mx=mr
            yield mx
    def BetaMatrix(self):
        m11 =self.matrix[0,0];  m12 =self.matrix[0,1]
        m21 =self.matrix[1,0];  m22 =self.matrix[1,1]
        n11 =self.matrix[2,2];  n12 =self.matrix[2,3]
        n21 =self.matrix[3,2];  n22 =self.matrix[3,3]
        m_beta = NP.array([
            [ m11*m11, -2.*m11*m12,           m12*m12,   0., 0., 0.],
            [-m11*m21,     m11*m22+m12*m21,  -m22*m12,   0., 0., 0.],
            [ m21*m21, -2.*m22*m21,           m22*m22,   0., 0., 0.],
            [ 0., 0., 0., n11*n11, -2.*n11*n12,           n12*n12],
            [ 0., 0., 0.,-n11*n21,     n11*n22+n12*n21,  -n22*n12],
            [ 0., 0., 0., n21*n21, -2.*n22*n21,           n22*n22]
            ])
        return m_beta
class I(_matrix):      ## unity matrix (an alias to _matrix class)
    def __init__(self,
    label='I',
    viseo=0.,
    beam=Beam.soll):
        super(I,self).__init__()
        self.label=label
        self.viseo=viseo
        self.beam=copy(beam)  # keep a local copy of the Beam instance (IMPORTANT!)
class D(I):            ## drift space nach Trace3D
    def __init__(self,
    length=0.,
    viseo=0.,
    label='D',
    beam=Beam.soll):    
        super(D,self).__init__(viseo=viseo,beam=beam)
        self.label=label
        self.length=length     ## hard edge length [m]
        self.matrix[0,1]=self.matrix[2,3]=self.length
        g=self.beam.gamma
        self.matrix[4,5]=length/(g*g)
    def shorten(self,l=0.):    # returns a new instance!
        return D(length=l,label=self.label,beam=self.beam,viseo=self.viseo)
    def update(self):          # returns a new instance!
        soll = Beam.soll
        return D(length=self.length,label=self.label,beam=soll,viseo=self.viseo)
class QF(D):           ## focusing quad nach Trace3D
    def __init__(self,
    k0=0.,
    length=0.,
    label='QF',
    beam=Beam.soll): 
        super(QF,self).__init__(length=length,label=label,beam=beam)
        self.k0=k0         ## Quad strength [m**-2]
        self.matrix=self._mx()
        self.viseo = +0.5
    def shorten(self,l=0.):
        return QF(k0=self.k0,length=l,label=self.label,beam=self.beam)
    def _mx(self):
        m=self.matrix
        g=self.beam.gamma
        rzz12=self.length/(g*g)
        kwurz=sqrt(self.k0)
        phi=self.length*kwurz
        ## focusing
        cf  =cos(phi)
        sf  =sin(phi)/kwurz
        cfp =-kwurz*sin(phi)
        sfp =cf
        ## defocusing
        cd  =cosh(phi)
        sd  =sinh(phi)/kwurz
        cdp =kwurz*sinh(phi)
        sdp =cd
        ## 6x6 matrix
        if (isinstance(self,QF)  and (isinstance(self,QD)==False)):
            m[0,0]=cf; m[0,1]=sf; m[1,0]=cfp; m[1,1]=sfp; m[2,2]=cd; m[2,3]=sd; m[3,2]=cdp; m[3,3]=sdp; m[4,5]=rzz12
        elif isinstance(self,QD):
            m[0,0]=cd; m[0,1]=sd; m[1,0]=cdp; m[1,1]=sdp; m[2,2]=cf; m[2,3]=sf; m[3,2]=cfp; m[3,3]=sfp; m[4,5]=rzz12
        else:
            raise RuntimeError('QF._mx: neither QF nor QD! should never happen!')
        return m
    def update(self):
        k0   =self.k0
        len  =self.length
        label=self.label
        beam =self.beam
        tki  =beam.tkin
        soll =Beam.soll
        tkf  =soll.tkin
        kf   =scalek0(k0,tki,tkf)  # scale quad strength
        # print('kf',kf)
        scaled=QF(k0=kf,length=len,label=label,beam=soll)
        return scaled        
class QD(QF):          ## defocusing quad nach Trace3D
    def __init__(self,
    k0=0.,
    length=0.,
    label='QD',
    beam=Beam.soll):
        super(QD,self).__init__(k0=k0,length=length,label=label,beam=beam)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0,length=l,label=self.label,beam=self.beam)
    def update(self):
        k0   =self.k0
        len  =self.length
        label=self.label
        beam =self.beam
        tki  =beam.tkin
        soll =Beam.soll
        tkf  =soll.tkin
        kf   =scalek0(k0,tki,tkf)
        scaled=QD(k0=kf,length=len,label=label,beam=soll)
        return scaled
class SD(D):           ## sector bending dipole in x-plane nach Trace3D
    def __init__(self,
    radius=0.,
    length=0.,
    label='SB',
    beam=Beam.soll):
        super(SD,self).__init__(length=length,label=label,beam=beam)
        self.radius = radius
        self.matrix=self._mx()
        self.viseo = 0.25
    def shorten(self,l=0.):
        return SD(radius=self.radius,length=l,label=self.label,beam=self.beam)
    def _mx(self):  # nach Trace3D
        m = self.matrix
        rho=self.radius
        k=1./rho
        phi=self.length/rho
        cx=cos(phi) ; sx=sin(phi)
        b=self.beam.beta
        ## x-plane
        m[0,0] = cx;     m[0,1] = sx/k;  m[0,5] = rho*(1.-cx)
        m[1,0] = -sx*k;  m[1,1] = cx;    m[1,5] = sx
        ## y-plane
        m[2,3] = self.length
        ## z-plane
        m[4,0] = -sx;   m[4,1] = -rho*(1.-cx);   m[4,5] = rho*sx-self.length*b*b
        return m
    def update(self):
        raise RuntimeWarning('SD.update(): not ready!')    
class RD(SD):          ## rectangular bending dipole in x-plane
    def __init__(self, 
    radius=0., 
    length=0., 
    label='RB',
    beam=Beam.soll):
        super(RD,self).__init__(radius=radius,length=length,label=label,beam=beam)
        wd = WD(self,label='',beam=beam)  # wedge myself...
        rd = wd * self * wd
        self.matrix= rd.matrix
    def shorten(self,l=0.):
        return RD(radius=self.radius,length=l,label=self.label,beam=self.beam)
    def update(self):
        raise RuntimeWarning('RD.update(): not ready!')    
class WD(D):           ## wedge of rectangular bending dipole in x-plane nach Trace3D
    def __init__(self,
    sector,
    label='WD',
    beam=Beam.soll):
        super(WD,self).__init__(label=label,beam=beam)
        m=self.matrix
        self.parent = sector
        self.radius = sector.radius
        self.psi = sector.length/self.radius
        rinv=1./self.radius
        psi=0.5*self.psi  ## Kantenwinkel
        ckp=rinv*tan(psi)
        ## 6x6 matrix
        m[1,0]=ckp
        m[3,2]=-ckp
    def shorten(self,l=0.):
        wd = WD(self.parent,label=self.label,beam=self.beam)
        m=wd.matrix
        wd.psi = l/wd.radius
        rinv=1./wd.radius
        psi=0.5*wd.psi  ## Kantenwinkel 
        ckp=rinv*tan(psi)
        ## 6x6 matrix
        m[1,0]=ckp
        m[3,2]=-ckp
        return wd
    def update(self):
        raise RuntimeWarning('WD.update(): not ready!')    
class CAV(D):          ## simple thin lens gap nach Dr.Tiede & T.Wrangler
    def __init__(self, 
    U0         =CONF['spalt_spannung'], 
    PhiSoll    =radians(CONF['soll_phase']), 
    fRF        =CONF['frequenz'], 
    label      ='RFG', 
    beam       =Beam(CONF['injection_energy']),
    gap        =CONF['spalt_laenge'],
    dWf=1.):
        super(CAV,self).__init__(label=label,beam=beam)
        self.u0     = U0                       # [MV] gap Voltage
        self.phis   = PhiSoll                  # [radians] soll phase
        self.freq   = fRF                      # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.beam.beta)              # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf      # T.Wrangler pp.221
        beami       = self.beam                               # Beam @ entrance
        beamf       = Beam(beami.tkin+self.deltaW)            # Beam @ exit
        dWavg       = (beamf.tkin - beami.tkin)*0.5+beami.tkin# average energy
        beam_avg    = Beam(dWavg)                             # Beam @ average energy
        b           = beam_avg.beta                           # beta @ average energy
        g           = beam_avg.gamma                          # gamma @ average energy
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx(self.tr,b,g)                   # transport matrix
        Beam.soll.incTK(self.deltaW)                          # beam accelerated
        self.viseo  = 0.25
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*self.freq*self.gap / (beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def _mx(self,tr,b,g):   # cavity nach Dr.Tiede pp.33 (todo: nach Trace3D)
        m=self.matrix
        e0 = self.beam.e0
        cyp = cxp = -pi*self.u0*tr*sin(self.phis)/(e0*self.lamb*g*g*g*b*b*b)  # T.Wrangler pp. 196
        # print(u"CAV: \u0394x'/x= ",cxp)
        # print("CAV: dx'/x= ",cxp)
        m[1,0]=cxp
        m[3,2]=cyp
        return m
    def shorten(self,l=0.):
        return self
    def update(self):
        return CAV(U0=self.u0,PhiSoll=self.phis,fRF=self.freq,label=self.label,beam=Beam.soll,gap=self.gap,dWf=self.dWf)
class RFG(D):          ## zero length RF gap nach Trace3D
    def __init__(self, 
    U0         =CONF['spalt_spannung'], 
    PhiSoll    =radians(CONF['soll_phase']), 
    fRF        =CONF['frequenz'], 
    label      ='RFG', 
    beam       =Beam.soll,
    gap        =CONF['spalt_laenge'],
    dWf=1.):
        super(RFG,self).__init__(label=label,beam=beam)
        self.u0     = U0*dWf                                  # [MV] gap Voltage
        self.phis   = PhiSoll                                 # [radians] soll phase
        self.freq   = fRF                                     # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.beam.beta)
        self.deltaW = self.u0*self.tr*cos(self.phis)          # Trace3D
        beami       = self.beam                               # Beam @ entrance
        beamf       = Beam(beami.tkin+self.deltaW)            # Beam @ exit
        dWavg       = (beamf.tkin - beami.tkin)*0.5+beami.tkin# average energy
        beam_avg    = Beam(dWavg)                             # Beam @ average energy
        b           = beam_avg.beta                           # beta @ average energy
        g           = beam_avg.gamma                          # gamma @ average energy
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx(self.tr,b,g,beami,beamf)       # transport matrix            
        Beam.soll.incTK(self.deltaW)
        # objprnt(Beam.soll,'soll')
        self.viseo  = 0.25      
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*self.freq*self.gap / (beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def _mx(self,tr,b,g,beami,beamf):   # RF gap nach Trace3D pp.17 (LA-UR-97-886)
        m=self.matrix
        e0 = self.beam.e0
        kz = 2.*pi*self.u0*tr*sin(self.phis)/(e0*b*b*self.lamb)  
        ky = kx = -0.5*kz/(g*g)
        bgi         = sqrt(beami.gamma*beami.gamma-1.)  # beta*gamma (i)
        bgf         = sqrt(beamf.gamma*beamf.gamma-1.)  # beta*gamma (f)
        bgiRbgf = bgi/bgf
        # 6x6
        m[1,0]=kx/bgf;    m[1,1]=bgiRbgf
        m[3,2]=ky/bgf;    m[3,3]=bgiRbgf
        m[5,4]=kz/bgf;    m[5,5]=bgiRbgf
        return m
    def shorten(self,l=0.):
        return self
    def update(self):
        return RFG(U0=self.u0,PhiSoll=self.phis,fRF=self.freq,label=self.label,beam=Beam.soll,gap=self.gap,dWf=self.dWf)
class _thin(_matrix):  ## the mother of all thin elements
    def __init__(self,beam=Beam.soll):
        self.beam = copy(beam)      ## keep a local copy of the Beam instance (important!)
    def step_through(self,anz=10):  ## stepping routine through the triplet (D,Kick,D)
        anz1 = int(anz/2)
        anz2 = int(anz-anz1*2)
        anz2 = int(anz1+anz2)
        for count,typ in enumerate(self.triplet):
            if count == 0:
                for i in range(anz1):
                    mx=typ.shorten(typ.length/anz1)
                    yield mx
            elif count == 1:
                mx=typ
                yield mx
            elif count == 2:
                for i in range(anz2):
                    mx=typ.shorten(typ.length/anz2)
                    yield mx
class QFth(_thin):     ## thin F-quad
    def __init__(self,
    k0=0.,
    length=0.,
    label='QFT',
    beam=Beam.soll):
        super(QFth,self).__init__(beam=beam)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,beam=self.beam,label=self.label,viseo=+0.5)
        df = di
        kick = _matrix()    ## 6x6 unit matrix
        m = kick.matrix     ## thin lens quad matrix
        if(self.k0l == 0.):
            m[1,0] = m[3,2] = 0.
        else:
            m[1,0] = -1./self.k0l
            m[3,2] = -m[1,0]
        lens = (di * kick) * df
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = +0.5
    def shorten(self,l=0.):
        raise RuntimeWarning('QFth.shorten(): not needed!')    
    def update(self):
        raise RuntimeWarning('QFth.update(): not ready!')    
class QDth(_thin):     ## thin D-quad
    def __init__(self,
    k0=0.,
    length=0.,
    label='QDT',
    beam=Beam.soll):
        super(QDth,self).__init__(beam=beam)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,beam=self.beam,label=self.label,viseo=-0.5)
        df = di
        kick = _matrix()    ## 6x6 unit matrix
        m = kick.matrix     ## thin lens quad matrix
        if(self.k0l == 0.):
            m[1,0] = m[3,2] = 0.
        else:
            m[1,0] = 1./self.k0l
            m[3,2] = -m[1,0]
        lens = (di * kick) * df
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = -0.5
    def shorten(self,l=0.):
        raise RuntimeWarning('QDth.shorten(): not needed!')    
    def update(self):
        raise RuntimeWarning('QDth.update(): not ready!')    
class RFC(_thin):      ## RF cavity as D*RFG*D
    def __init__(self,
    U0=CONF['spalt_spannung'],
    PhiSoll=CONF['soll_phase'],
    fRF=CONF['frequenz'],
    label='RFC',
    beam=Beam.soll,
    gap=CONF['spalt_laenge'],
    length=0.,
    dWf=1.):
        super(RFC,self).__init__(beam=beam)
        self.u0     = U0*dWf
        self.phis   = PhiSoll
        self.freq   = fRF
        self.label  = label
        self.gap    = gap
        self.length = length
        self.dWf    = dWf
        di   = D(length=0.5*length,label='D(i)',beam=beam)
        kick = RFG(U0=self.u0,PhiSoll=self.phis,fRF=self.freq,label=self.label,beam=self.beam,gap=self.gap,dWf=self.dWf)  ## Trace3D RF gap
        self.tr     = kick.tr
        # objprnt(kick)
        df   = D(length=0.5*length,label='D(f)',beam=Beam.soll)   # energy update here
        lens = (di * kick) * df
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
    def shorten(self,l=0.):
        raise RuntimeWarning('RFC.shorten(): not needed!')    
    def update(self):
        return RFC(U0=self.u0,PhiSoll=self.phis,fRF=self.freq,label=self.label,beam=Beam.soll,gap=self.gap,length=self.length,dWf=self.dWf)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
class Test(_matrix):
    def __init__(self,a,b,c,d,e,f,label='test'):
        super(Test,self).__init__()
        self.matrix=NP.array([[a,b,0.,0.,0.,0.],
                              [c,d,0.,0.,0.,0.],
                              [0.,0.,a,b,0.,0.],
                              [0.,0.,d,e,0.,0.],
                              [0.,0.,0.,0.,a,b],
                              [0.,0.,0.,0.,e,f]])
        self.label=label
def k0test(gradient=0.,beta=0.,energy=0.):   ## helper function for tests
    """
        quad strength as function of energy and gradient
        gradient in [Tesla/m]
        energy in [Gev]
        beta [v/c]
    """
    if (gradient != 0. and energy != 0. and beta !=0.):
        return 0.2998*gradient/(beta*energy)             
    else:
        raise RuntimeError('zero gradient or energy or beta in quad strength!')
def test0():
    print('trivial test 0 ...')
    a=Test(1,2,3,4,5,6,label='a')
    a.out()
    b=Test(1,1,1,1,1,1,label='b')
    b.out()
    (a*b).out()
    (b*a).out()
    print('--------------- EOF test0 --------------------')
def test1():
    print('trivial test 1 ...')
    i1=_matrix()
    i2=i1*i1
    i1.out()
    i2.out()
    print('--------------- EOF test1 --------------------')
def test2():
    print('trivial test 2 ...')
    i1=_matrix()
    d1=D(10.,'D1')
    d1.out()
    (d1*d1).out()
    (d1*d1*d1).out()
    (d1*i1*d1).out()
    d2=D(90,'D2')
    d2.out()
    (d1*d2).out()
    d3=D(90.,label='')
    (d2*d3).out()
    print('--------------- EOF test2 --------------------')
def test3():
    print('test product of _matrix class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.)
    qf.out()
    ## test product of _matrix class
    qd=QD(k0=k,length=1.)
    qd.out()
    (qf*qd).out()
    print('--------------- EOF test3 --------------------')
def test4():
    print('test shortening of elements ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    ## elements
    d10=D(10.,'d10')
    d10.out()
    (d10.shorten(1.e-2)).out()
    
    qf=QF(k0=k,length=1.)
    qf05=qf.shorten(0.2)
    (qf05*qf05*qf05*qf05*qf05).out()
    qf.out()
    
    qd=QD(k0=k,length=1.)
    qd05=qd.shorten(0.2)
    (qd05*qd05*qd05*qd05*qd05).out()
    qd.out()
    
    sd=SD(radius=10.,length=2.)
    sd05=sd.shorten(0.4)
    (sd05*sd05*sd05*sd05*sd05).out()
    sd.out()
    print('--------------- EOF test4 --------------------')
def test5():
    print("K.Wille's Beispiel auf pp. 112-113")
    kqf=  wille()['k_quad_f']
    lqf=  wille()['length_quad_f']
    kqd=  wille()['k_quad_d']
    lqd=  wille()['length_quad_d']
    rhob= wille()['beding_radius']
    lb=   wille()['dipole_length']
    ld=   wille()['drift_length']
    ## elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(mb)
    md=D(ld)
    ## test matrix multiplication
    mz=I()
    mz=mz *mqf
    mz=mz *md
    mz=mz *mw
    mz=mz *mb
    mz=mz *mw
    mz=mz *md
    mz=mz *mqd
    mz=mz *md
    mz=mz *mw
    mz=mz *mb
    mz=mz *mw
    mz=mz *md
    mz=mz *mqf
    mz.out()
    print('--------------- EOF test5 --------------------')
def test6():
    print('test step_through elements ...')
    kqf=  wille()['k_quad_f']
    lqf=  wille()['length_quad_f']
    kqd=  wille()['k_quad_d']
    lqd=  wille()['length_quad_d']
    rhob= wille()['beding_radius']
    lb=   wille()['dipole_length']
    ld=   wille()['drift_length']

    ## elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(mb)
    md=D(ld)
    
    ## test step_through elements ...
    list=[mqf,mqd,mb,mw,md]
    # list=[mqf]
    for m_anfang in list:
        m_end=I()
        print('======================================')
        for count,mi in enumerate(m_anfang.step_through()):
            print('step ',count+1,end='  ')
            mi.out()
            m_end=m_end*mi
        m_end.out()
        m_anfang.out()
    print('--------------- EOF test6 --------------------')
def test7():
    print('======================================')
    print('test Rechteckmagnet...')
    rhob= wille()['beding_radius'] 
    lb=   wille()['dipole_length']
    mb=SD(radius=rhob,length=lb,label='B')
    mw=WD(mb,label='W')
    mr=mw*mb*mw
    mw.out()
    mb.out()
    mr.out()
    mr=RD(radius=rhob,length=lb,label='R')
    mr.out()
    print('--------------- EOF test7 --------------------')
def test8():
    print('test cavity...')
    objprnt(Beam.soll,'soll')
    cav=CAV()
    objprnt(cav,'CAV')
    objprnt(Beam.soll,'soll')
    rfg=RFG()
    objprnt(rfg,'RFG')
    objprnt(Beam.soll,'soll')
    print('--------------- EOF test8 --------------------')
def test9():
    print('\ntest: quad k-faktor and quad scaling')
    grad=CONF['quad_gradient']   # [T/m] gradient
    tk=CONF['injection_energy']    # [MeV]  kin. energy
    kq=k0(gradient=grad,tkin=tk)    # quad strength [1/m**2]
    len=0.4                         # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('\nproton {:.3f}[MeV] ~ beta {:.3f} in quadrupole:'.format(tk,Beam(tk).beta))
    print('k [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))
    
    grad=dBdz(kq,tk) # quad gradient from k and tkinetic
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    
    mqf=QF(kq,len)
    mqd=QD(kq,len)
    cavity=CAV(
        U0=CONF['spalt_spannung'],
        PhiSoll=radians(CONF['soll_phase']),
        fRF=CONF['frequenz'],
        label='gap')
    print('========================')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    Beam.soll = Beam(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        Beam.soll.incTK(dt)
        mqf.update().out()
    print('========================')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    Beam.soll = Beam(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        Beam.soll.incTK(dt)
        mqd.update().out()
    print('========================')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    Beam.soll = Beam(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        Beam.soll.incTK(dt)
        cavity.update().out()
    print('--------------- EOF test9 --------------------')
def test10():
    print('\ntest: Beam class')
    dictprnt(CONF,text='setup.CONF')

    # Beam class
    print()
    Beam(0.).out()
    Beam(50.).out()
    Beam(200.).out()
    Beam(1.e6).out()
    Beam(1.e9).out()
    
    beam = Beam.soll
    beam.out()
    beam.incTK(150.)
    beam.out()
    print('--------------- EOF test10 --------------------')
def test11():
    print('\ntest thin lenses:')
    print('----------------- product matrix ---------')
    k0=1.
    length = 2.
    qf=QFth(k0=k0,length=length)
    objprnt(qf,'QFthin')
    qd=QDth(k0=k0,length=length)
    objprnt(qd,'QDthin')
    qf.out()
    qd.out()
    print('---------------- step through ---------------')
    for elm in qf.step_through(6):
        elm.out()
    for elm in qd.step_through(7):
        elm.out()
    print('------ RF cavity test & step through --------')
    objprnt(Beam.soll,'beam(i)')
    rf=RFC(length=length)
    rf.out()
    objprnt(Beam.soll,'beam(f)')
    objprnt(rf,'RFC cavity')
    for elm in rf.step_through():
        elm.out()
    print('--------------- EOF test11 --------------------')
if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()
    test11()
