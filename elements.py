# -*- coding: utf-8 -*-
from setup import wille, Phys
import numpy as NP 
from math import sqrt, sinh, cosh, sin, cos, fabs, tan, floor, modf, pi

class Beam():   ## relativistic particles
    def __init__(self,tkin=0.):
        self._set_self(tkin)
    def _set_self(self,tkin):
        self.tkin = tkin                     # proton kinetic energy [MeV]
        self.e0   = Phys['proton_mass']      # proton rest mass [MeV/c**2]
        self.e    = self.e0+self.tkin        # proton total energy [MeV]
        self.gamma= self.e/self.e0
        self.beta = sqrt(1.-1./(self.gamma*self.gamma))
        self.v    = self.beta*Phys['lichtgeschwindigkeit']
        self.name = 'proton'
    def incTK(self,deltaTK):
        self._set_self(self.tkin+deltaTK)
    def out(self):
        print('{:s}:  T-kin[MeV]={:.3f} gamma {:.3f} beta {:.3f} velocity[m/s] {:.6g} E[MeV] {:.3f} '
            .format(self.name,self.tkin,self.gamma,self.beta,self.v,self.e))        
Beam.soll = Beam(Phys['kinetic_energy']) # the synchronous reference particle  (class member!)
class Matrix(object):  # 6x6 matrices
    _dim = 6   # 6x6 matrices
    def __init__(self):
        self.matrix=NP.eye(Matrix._dim)    ## 6x6 unit matrix
        self.label=''
        self.length=0.         ## default zero length!
        self.slice_min = 0.01  ## minimal slice length
        self.viseo = 0.
    def out(self):
        print(self.label)
        print(self.matrix)
    def __mul__(self,other):
        product=NP.einsum('ij,jk',self.matrix,other.matrix)
        res=Matrix()
        if (self.label == ''):
            res.label=other.label
        else:
            res.label=self.label+'*'+other.label
        res.length=self.length+other.length
        res.matrix=product
        return res
    def reverse(self):
        res=Matrix()
        for i in range(Matrix._dim):
            for k in range(Matrix._dim):
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
        return Matrix()
    def step_through(self,anz=10):
        """
        Step through an element(The central nontrivial function).
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
class I(Matrix):## unity Matrix (an alias to Matrix class)
    def __init__(self,label='I',viseo=0.,beam=Beam(Phys['kinetic_energy'])):
        super(I,self).__init__()
        self.label=label
        self.viseo=viseo
        self.beam=beam
class Test(Matrix):
    def __init__(self,a,b,c,d,e,f,label='test'):
        super(Test,self).__init__()
        self.matrix=NP.array([[a,b,0.,0.,0.,0.],
                              [c,d,0.,0.,0.,0.],
                              [0.,0.,a,b,0.,0.],
                              [0.,0.,d,e,0.,0.],
                              [0.,0.,0.,0.,a,b],
                              [0.,0.,0.,0.,e,f]])
        self.label=label
class D(Matrix):## drift space nach Trace3D
    def __init__(self,length=0.,label='D',beam=Beam(Phys['kinetic_energy'])):    
        super(D,self).__init__()
        self.label=label
        self.beam=beam
        self.length=length     ## hard edge length [m]
        self.matrix[0,1]=self.matrix[2,3]=self.length
        g=self.beam.gamma
        self.matrix[4,5]=self.length/(g*g)
    def shorten(self,l=0.):
        return D(length=l,label=self.label,beam=self.beam)
class QF(D):    ## focusing quad nach Trace3D
    def __init__(self,k0=0.,length=0.,label='QF',beam=Beam(Phys['kinetic_energy'])):    
        super(QF,self).__init__(length=length,label=label,beam=beam)
        self.k0=k0         ## energy independent Quad strength [m**-2]
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
        elif isinstance(self,QD) :
             m[0,0]=cd; m[0,1]=sd; m[1,0]=cdp; m[1,1]=sdp; m[2,2]=cf; m[2,3]=sf; m[3,2]=cfp; m[3,3]=sfp; m[4,5]=rzz12
        else:
            raise RuntimeError('QF._mx: neither QF nor QD! should never happen!')
        return m
    def scale(self,deltaTk=0.):
        # for k,v in self.__dict__.items():
            # print(k.rjust(30),':',v)
        # for k,v in self.beam.__dict__.items():
            # print(k.rjust(30),':',v)
        k0   =self.k0
        len  =self.length
        label=self.label
        tki = self.beam.tkin
        tkf = tki+deltaTk
        kf=scalek0(k0,tki,tkf)
        # print('kf',kf)
        quad_scaled=QF(k0=kf,length=len,label=label,beam=self.beam)
        return quad_scaled
class QD(QF):   ## defocusing quad nach Trace3D
    def __init__(self,k0=0.,length=0.,label='QD',beam=Beam(Phys['kinetic_energy'])):
        super(QD,self).__init__(k0=k0,length=length,label=label,beam=beam)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0,length=l,label=self.label,beam=self.beam)
    def scale(self,deltaTk=0.):
        k0   =self.k0
        len  =self.length
        label=self.label
        tki = self.beam.tkin
        tkf = tki+deltaTk
        kf=scalek0(k0,tki,tkf)
        quad_scaled=QD(k0=kf,length=len,label=label,beam=self.beam)
        return quad_scaled
class SD(D):    ## sector bending dipole in x-plane nach Trace3D
    def __init__(self,radius=0.,length=0.,label='SB',beam=Beam(Phys['kinetic_energy'])):
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
class RD(SD):   ## rectangular bending dipole in x-plane
    def __init__(self, radius=0., length=0., label='RB',beam=Beam(Phys['kinetic_energy'])):
        super(RD,self).__init__(radius=radius,length=length,label=label,beam=beam)
        wd = WD(self,label='',beam=beam)  # wedge myself...
        rd = wd * self * wd
        self.matrix= rd.matrix
    def shorten(self,l=0.):
        return RD(radius=self.radius,length=l,label=self.label,beam=self.beam)
class WD(D):    ## wedge of rectangular bending dipole in x-plane nach Trace3D
    def __init__(self,sector,label='WD',beam=Beam(Phys['kinetic_energy'])):
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
class CAV(D):   ## simple thin lens gap nach Dr.Tiede & T.Wrangler
    def __init__(self,U0=10.,PhiSoll=-0.25*pi,fRF=800.,label='CAV',beam=Beam(Phys['kinetic_energy']),dWf=1.):
        super(CAV,self).__init__(label=label,beam=beam)
        self.u0     = U0       # [MV] gap Voltage
        self.phis   = PhiSoll  # [radians] soll phase
        self.freq   = fRF      # [MHz]  RF frequenz
        self.lamb   = 1.e-6*Phys['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.beam.beta) # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf # T.Wrangler pp.221
        beami       = self.beam                      # Beam @ entrance
        beamf       = Beam(beami.tkin+self.deltaW)   # Beam @ exit
        dWavg       = (beamf.tkin - beami.tkin)*0.5+beami.tkin  # average energy
        beam_avg    = Beam(dWavg)       # Beam @ average energy
        b           = beam_avg.beta     # beta @ average energy
        g           = beam_avg.gamma    # gamma @ average energy
        self.Ks     = 2.*pi/(self.lamb*g*b)  # T.Wrangler pp.196
        self.matrix = self._mx(self.tr,b,g)  # transport matrix
        self.beam   = beamf             # Beam @ exit
        sollParticle= beamf
        Beam.soll.incTK(self.deltaW)
        self.viseo  = 0.25
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        gap_len = Phys['spalt_laenge']
        teta = 2.*pi*1.e6*self.freq*gap_len / (beta*Phys['lichtgeschwindigkeit'])
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
    def scale(self,deltaTk=0.):
        u0=self.u0
        phiSoll=self.phis
        fRF=self.freq
        label=self.label
        tki =self.beam.tkin
        tkf =tki+deltaTk
        cav = CAV(U0=u0,PhiSoll=phiSoll,fRF=fRF,beam=Beam(tkf),label=label)
        return cav
class RFG(D):   ## thin lens RF gap nach Trace3D
    def __init__(self, U0=10., PhiSoll=-0.25*pi, fRF=800., label='RFG', beam=Beam(Phys['kinetic_energy']),dWf=1.):
        super(RFG,self).__init__(label=label,beam=beam)
        self.u0     = U0       # [MV] gap Voltage
        self.phis   = PhiSoll  # [radians] soll phase
        self.freq   = fRF      # [MHz]  RF frequenz
        self.lamb   = 1.e-6*Phys['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.beam.beta)
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf # Trace3D
        beami       = self.beam                      # Beam @ entrance
        beamf       = Beam(beami.tkin+self.deltaW)   # Beam @ exit
        dWavg       = (beamf.tkin - beami.tkin)*0.5+beami.tkin  # average energy
        beam_avg    = Beam(dWavg)       # Beam @ average energy
        b           = beam_avg.beta     # beta @ average energy
        g           = beam_avg.gamma    # gamma @ average energy
        self.matrix = self._mx(self.tr,b,g,beami,beamf)   # transport matrix            
        self.beam   = beamf             # Beam @ exit
        Beam.soll.incTK(self.deltaW)
        # dictp(Beam.soll,'soll')
        self.viseo  = 0.25      
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        gap_len = Phys['spalt_laenge']
        teta = 2.*pi*1.e6*self.freq*gap_len / (beta*Phys['lichtgeschwindigkeit'])
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
    def scale(self,deltaTk=0.):
        u0=self.u0
        phiSoll=self.phis
        fRF=self.freq
        label=self.label
        tki = self.beam.tkin
        tkf = tki+deltaTk
        rfg = RFG(U0=u0,PhiSoll=phiSoll,fRF=fRF,beam=Beam(tkf),label=label)
        return rfg
def k0(gradient=0.,tkin=0.):       ## quad strength from B-field gradient & kin. energy
    """
    quad strength as function of kin. energy and gradient
    gradient: in [Tesla/m]
    tkin: kinetic energy in [MeV]
    """
    if (tkin >= 0.):
        prot=Beam(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*Phys['lichtgeschwindigkeit']
        kres = factor*gradient/(beta*gamma*e0) 
        # print('k0= ',kres)
        return kres
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
def scalek0(k0=0.,tki=0.,tkf=0.):  ## scale quad  strength with kin. energy
    """
    scale k0 for increase of kin. energy from
    tki to tkf
    """
    pi  =Beam(tki)
    bi  =pi.beta
    gi  =pi.gamma
    pf  =Beam(tkf)
    bf  =pf.beta
    gf  =pf.gamma
    kf= k0 * (bi * gi) / (bf * gf)
    return kf
def dBdz(k0=0.,tkin=0.):           ## B-field gradient from quad strength & kin. energy
    """
    calculate quad gradient for given quad strength k0
    and given kin. energy tkin
    """
    if (tkin >= 0.):
        prot=Beam(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*Phys['lichtgeschwindigkeit']
        return k0*(beta*gamma*e0)/factor           
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
####################################################################
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
def dictp(what,text='========',filter={}):   ## helper
        print('========= '+text+' =================')
        for k,v in what.__dict__.items():
            if k in filter:
                continue
            print(k.rjust(30),':',v)
def test0():
    print('trivial test 0 ...')
    a=Test(1,2,3,4,5,6,label='a')
    a.out()
    b=Test(1,1,1,1,1,1,label='b')
    b.out()
    (a*b).out()
    (b*a).out()
def test1():
    print('trivial test 1 ...')
    i1=Matrix()
    i2=i1*i1
    i1.out()
    i2.out()
def test2():
    print('trivial test 2 ...')
    i1=Matrix()
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
def test3():
    print('test product of Matrix class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.)
    qf.out()
    ## test product of Matrix class
    qd=QD(k0=k,length=1.)
    qd.out()
    (qf*qd).out()
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
def test8():
    print('test cavity...')
    dictp(Beam.soll,'soll')
    cav=CAV()
    dictp(cav,'CAV')
    dictp(Beam.soll,'soll')
    rfg=RFG()
    dictp(rfg,'RFG')
    dictp(Beam.soll,'soll')
def test9():
    print('\ntest: quad k-faktor and quad scaling')
    grad=Phys['quad_gradient']   # [T/m] gradient
    tk=Phys['kinetic_energy']    # [MeV]  kin. energy
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
        U0=Phys['spalt_spannung'],
        PhiSoll=Phys['soll_phase']*Phys['radians'],
        fRF=Phys['frequenz'],
        beam=Beam(tk),
        label='gap')
    tki=tk
    for dt in [10.,50.,150.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        mqf.scale(deltaTk=dt).out()
        mqd.scale(deltaTk=dt).out()
        cavity.scale(deltaTk=dt).out()
def test10():
    print('\ntest: Beam class')
    for key,value in Phys.items():
        print('{}=  {:.4g}'.format(key.rjust(20),value))
    # Beam class
    print()
    Beam(0.).out()
    Beam(50.).out()
    Beam(200.).out()
    Beam(1.e6).out()
    Beam(1.e9).out()
    
    beam = sollParticle
    beam.out()
    beam.incTK(150.)
    beam.out()
if __name__ == '__main__':
    # test0()
    # test1()
    test2()
    # test3()
    # test4()
    # test5()
    # test6()
    # test7()
    # test8()
    # test9()
    # test10()
