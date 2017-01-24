#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
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
from math import sqrt,sinh,cosh,sin,cos,fabs,tan,floor,modf,pi,radians,degrees
from copy import copy
import numpy as NP

from setutils import wille,CONF,dictprnt,objprnt,Particle
from setutils import DEBUG,k0,dBdz,scalek0,printv

MDIM=10        # dimension of matrices

# x        x'        y        y'       z         z'       E        dE        s        l       dispersion
XKOO = 0;XPKOO = 1;YKOO = 2;YPKOO = 3;ZKOO = 4;ZPKOO = 5;EKOO = 6;DEKOO = 7;SKOO = 8;LKOO = 9;#DISP = 10

NP.set_printoptions(linewidth=132,formatter={'float':'{:>8.5g}'.format})  #pretty printing

class _matrix(object): ## the mother of all matrices
    # MDIMxMDIM matrices used here
    def __init__(self):
        self.matrix=NP.eye(MDIM)    ## MDIMxMDIM unit matrix
        self.label=''
        self.length=0.         ## default zero length!
        self.slice_min = 0.01  ## minimal slice length
        self.viseo = 0.
    def string(self):
        s='{}\n'.format(self.label)
        for i in range(MDIM):
            for j in range(MDIM):
                s+='{:8.4g} '.format(self.matrix[i,j])
            s+='\n'
        return s
    def __mul__(self,other):
#         DEBUG('A*B >>',self.label,'*',other.label)
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
        raise RuntimeError('_matrix:reverse not released yet!')
        res=_matrix()
        for i in range(MDIM):
            for k in range(MDIM):
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
    def shorten(self,length=0.):    # virtual function to be implemented by child classes
        raise RuntimeError('_matrix.shorten(): empty virtual member function called!')
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
            (rest,fanz) = modf(self.length/step)
            anz = int(fanz)
            rest = self.length * rest
            mx = self.shorten(step)
            if fabs(rest) > 1.e-9:
                mr = self.shorten(rest)
            else:
                mr=I(label=self.label,viseo=self.viseo)
        # DEBUG('label={} fanz={} anz={} step={} rest={}'.format(self.label,fanz,anz,step,rest))
        for i in range(anz+1):
            if i == anz:
                mx=mr
            yield mx
    def betaMatrix(self):
#         m11 =self.matrix[0,0];  m12 =self.matrix[0,1]
#         m21 =self.matrix[1,0];  m22 =self.matrix[1,1]
#         n11 =self.matrix[2,2];  n12 =self.matrix[2,3]
#         n21 =self.matrix[3,2];  n22 =self.matrix[3,3]
        m11 =self.matrix[XKOO,XKOO];   m12 =self.matrix[XKOO,XPKOO]
        m21 =self.matrix[XPKOO,XKOO];  m22 =self.matrix[XPKOO,XPKOO]
        n11 =self.matrix[YKOO,YKOO];   n12 =self.matrix[YKOO,YPKOO]
        n21 =self.matrix[YPKOO,YKOO];  n22 =self.matrix[YPKOO,YPKOO]
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
    def __init__(self, label='I', viseo=0., particle=Particle.soll):
        super(I,self).__init__()
        self.label=label
        self.viseo=viseo
        self.particle=copy(particle)  # keep a local copy of the particle instance (IMPORTANT!)

class MRK(I):          ## a marker
    def __init__(self, label='MRK', particle=Particle.soll):
        super(MRK,self).__init__(label=label, particle=particle)
    def shorten(self,l=0):
        return self
    def adapt_for_energy(self,tkin):
        self.__init__(label=self.label, particle=Particle(tkin=tkin))
        return self

class D(I):            ## drift space nach Trace3D
    def __init__(self, length=0., viseo=0., label='D', particle=Particle.soll):
        super(D,self).__init__(viseo=viseo,particle=particle)
        self.label=label
        self.length=length     ## hard edge length [m]
        g=self.particle.gamma
        self.matrix[XKOO,XPKOO] = self.matrix[YKOO,YPKOO] =self.length
        self.matrix[ZKOO,ZPKOO] = self.length/(g*g)
        self.matrix[SKOO,LKOO]  = self.length     #delta-s
    def shorten(self,l=0.):    # returns a new instance!
        return D(length=l,label=self.label,particle=self.particle,viseo=self.viseo)
    def adapt_for_energy(self,tkin):
        self.__init__(length=self.length, viseo=self.viseo, label=self.label, particle=Particle(tkin=tkin))
        return self

class QF(D):           ## focusing quad nach Trace3D
    def __init__(self, k0=0., length=0., label='QF', particle=Particle.soll):
        super(QF,self).__init__(length=length,label=label,particle=particle)
        self.k0=k0         ## Quad strength [m**-2]
        self.matrix=self._mx()
        self.viseo = +0.5
    def shorten(self,l=0.):
        return QF(k0=self.k0,length=l,label=self.label,particle=self.particle)
    def _mx(self):
        m=self.matrix
        g=self.particle.gamma
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
        #       0    0           0    1            1     0            1     1            2     2
            m[XKOO,XKOO]=cf; m[XKOO,XPKOO]=sf; m[XPKOO,XKOO]=cfp; m[XPKOO,XPKOO]=sfp; m[YKOO,YKOO]=cd
        #       2    3            3    2             3     3            4     5
            m[YKOO,YPKOO]=sd; m[YPKOO,YKOO]=cdp; m[YPKOO,YPKOO]=sdp; m[ZKOO,ZPKOO]=rzz12
        elif isinstance(self,QD):
            m[XKOO,XKOO]=cd; m[XKOO,XPKOO]=sd; m[XPKOO,XKOO]=cdp; m[XPKOO,XPKOO]=sdp; m[YKOO,YKOO]=cf; m[YKOO,YPKOO]=sf; m[YPKOO,YKOO]=cfp; m[YPKOO,YPKOO]=sfp; m[ZKOO,ZPKOO]=rzz12
        else:
            raise RuntimeError('QF._mx: neither QF nor QD! should never happen!')
        return m
    def adapt_for_energy(self,tkin):
        kf = scalek0(self.k0,self.particle.tkin,tkin)
        self.__init__(k0=kf, length=self.length, label=self.label, particle=Particle(tkin=tkin))
        return self

class QD(QF):          ## defocusing quad nach Trace3D
    def __init__(self, k0=0., length=0., label='QD', particle=Particle.soll):
        super(QD,self).__init__(k0=k0,length=length,label=label,particle=particle)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0,length=l,label=self.label,particle=self.particle)

class SD(D):           ## sector bending dipole in x-plane nach Trace3D
#
# ACHTUNG Matrix nicht vollstaedig; muss ueberprueft werden!
    def __init__(self,
                radius=0.,
                length=0.,
                label='SB',
                particle=Particle.soll):
        raise RuntimeError('SD:not implemented!')
        super(SD,self).__init__(length=length,label=label,particle=particle)
        self.radius = radius
        self.matrix=self._mx()
        self.viseo = 0.25
    def shorten(self,l=0.):
        return SD(radius=self.radius,length=l,label=self.label,particle=self.particle)
    def _mx(self):  # nach Trace3D
        m = self.matrix
        rho=self.radius
        k=1./rho
        phi=self.length/rho
        cx=cos(phi) ; sx=sin(phi)
        b=self.particle.beta
        ## x-plane
#         m[0,0] = cx;          m[0,1] = sx/k;           m[0,5] = rho*(1.-cx)
#         m[1,0] = -sx*k;       m[1,1] = cx;             m[1,5] = sx
        m[XKOO,XKOO]  = cx;     m[XKOO,XPKOO]   = sx/k;  m[XKOO,ZPKOO]  = rho*(1.-cx)
        m[XPKOO,XKOO] = -sx*k;  m[XPKOO,XPKOO]  = cx;    m[XPKOO,ZPKOO] = sx
        ## y-plane
#         m[2,3] = self.length
        m[YKOO,YPKOO] = self.length
        ## z-plane
#         m[4,0] = -sx;       m[4,1] = -rho*(1.-cx);          m[4,5] = rho*sx-self.length*b*b
        m[ZKOO,XKOO] = -sx;   m[ZKOO,XPKOO] = -rho*(1.-cx);   m[ZKOO,ZPKOO] = rho*sx-self.length*b*b
        return m

class RD(SD):          ## rectangular bending dipole in x-plane
# ACHTUNG Matrix nicht vollstaedig; muss ueberprueft werden!
    def __init__(self,
                radius=0.,
                length=0.,
                label='RB',
                particle=Particle.soll):
        raise RuntimeError('RD:not implemented!')
        super(RD,self).__init__(radius=radius,length=length,label=label,particle=particle)
        wd = WD(self,label='',particle=particle)  # wedge myself...
        rd = wd * self * wd
        self.matrix= rd.matrix
    def shorten(self,l=0.):
        return RD(radius=self.radius,length=l,label=self.label,particle=self.particle)

class WD(D):           ## wedge of rectangular bending dipole in x-plane nach Trace3D
# ACHTUNG Matrix nicht vollstaedig; muss ueberprueft werden!
    def __init__(self,
                sector,
                label='WD',
                particle=Particle.soll):
        raise RuntimeError('WD:not implemented!')
        super(WD,self).__init__(label=label,particle=particle)
        m=self.matrix
        self.parent = sector
        self.radius = sector.radius
        self.psi = sector.length/self.radius
        rinv=1./self.radius
        psi=0.5*self.psi  ## Kantenwinkel
        ckp=rinv*tan(psi)
        ## 6x6 matrix
#         m[1,0]=ckp
#         m[3,2]=-ckp
        m[XPKOO,XKOO]=ckp
        m[YPKOO,YKOO]=-ckp
    def shorten(self,l=0.):
        wd = WD(self.parent,label=self.label,particle=self.particle)
        m=wd.matrix
        wd.psi = l/wd.radius
        rinv=1./wd.radius
        psi=0.5*wd.psi  ## Kantenwinkel
        ckp=rinv*tan(psi)
        ## 6x6 matrix
        m[XPKOO,XKOO]=ckp
        m[YPKOO,YKOO]=-ckp
        return wd

class CAV(D):          ## simple thin lens gap nach Dr.Tiede & T.Wrangler
    def __init__(self,
                        U0         =CONF['spalt_spannung'],
                        PhiSoll    =radians(CONF['soll_phase']),
                        fRF        =CONF['frequenz'],
                        label      ='RFG',
                        particle   =Particle(CONF['injection_energy']),
                        gap        =CONF['spalt_laenge'],
                        dWf=1.):
        super(CAV,self).__init__(label=label,particle=particle)
        self.u0     = U0                       # [MV] gap Voltage
        self.phis   = PhiSoll                  # [radians] soll phase
        self.freq   = fRF                      # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.particle.beta)          # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf      # T.Wrangler pp.221
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        part_center = Particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap center
        g           = part_center.gamma                       # gamma @ gap center
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx(self.tr,b,g)                   # transport matrix
        self.viseo  = 0.25
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*self.freq*self.gap / (beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def _mx(self,tr,b,g):   # cavity nach Dr.Tiede pp.33 (todo: nach Trace3D)
        m=self.matrix
        e0 = self.particle.e0
        cyp = cxp = -pi*self.u0*tr*sin(self.phis)/(e0*self.lamb*g*g*g*b*b*b)  # T.Wrangler pp. 196
#         m[1,0]=cxp
#         m[3,2]=cyp
#         m[6,7]=self.deltaW      #energy kick = acceleration
        m[XPKOO,XKOO]=cxp
        m[YPKOO,YKOO]=cyp
        m[EKOO,DEKOO]=self.deltaW      #energy kick = acceleration
        return m
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        return self

class RFG(D):          ## zero length RF gap nach Trace3D
    def __init__(self,
                    U0         =CONF['spalt_spannung'],
                    PhiSoll    =radians(CONF['soll_phase']),
                    fRF        =CONF['frequenz'],
                    label      ='RFG',
                    particle   =Particle.soll,
                    gap        =CONF['spalt_laenge'],
                    dWf=1.):
        super(RFG,self).__init__(label=label,particle=particle)
        self.viseo  = 0.25
        self.u0     = U0*dWf                                  # [MV] gap Voltage
        self.phis   = PhiSoll                                 # [radians] soll phase
        self.freq   = fRF                                     # [Hz]  RF frequenz
        self.label  = label
        self.gap    = gap
        self.dWf    = dWf
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._TrTF(self.particle.beta)
        self.deltaW = self.u0*self.tr*cos(self.phis)          # Trace3D
#         DEBUG('\n',self.particle.string())
#         DEBUG('RFG U0,phis,tr >>',self.u0,degrees(self.phis),self.tr)
#         DEBUG('RFG deltaW >>',self.deltaW)
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        part_center = Particle(tk_center)                     # particle @ gap center
        particlei   = self.particle                           # particle @ entrance
        particlef   = Particle(particlei.tkin+self.deltaW)    # particle @ exit
        b           = part_center.beta                        # beta @ gap cemter
        g           = part_center.gamma                       # gamma @ gap center
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx(self.tr,b,g,particlei,particlef)       # transport matrix
    def _TrTF(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = pi*self.freq*self.gap / CONF['lichtgeschwindigkeit']
#         DEBUG('teta , beta>>',teta,beta)
        teta = teta/beta
        ttf = sin(teta)/teta
        return ttf
    def _mx(self,tr,b,g,particlei,particlef):   # RF gap nach Trace3D pp.17 (LA-UR-97-886)
        m=self.matrix
        e0 = self.particle.e0
        kz = 2.*pi*self.u0*tr*sin(self.phis)/(e0*b*b*self.lamb)
        ky = kx = -0.5*kz/(g*g)
        bgi         = sqrt(particlei.gamma*particlei.gamma-1.)  # beta*gamma (i)
        bgf         = sqrt(particlef.gamma*particlef.gamma-1.)  # beta*gamma (f)
        bgiRbgf = bgi/bgf
        # 6x6
#         m[1,0]=kx/bgf;    m[1,1]=bgiRbgf
#         m[3,2]=ky/bgf;    m[3,3]=bgiRbgf
#         m[5,4]=kz/bgf;    m[5,5]=bgiRbgf
        m[XPKOO,XKOO]=kx/bgf;    m[XPKOO,XPKOO]=bgiRbgf
        m[YPKOO,YKOO]=ky/bgf;    m[YPKOO,YPKOO]=bgiRbgf
        m[ZPKOO,ZKOO]=kz/bgf;    m[ZPKOO,ZPKOO]=bgiRbgf
        #energy kick at acceleration
#         m[6,7]=self.deltaW
        m[EKOO,DEKOO]=self.deltaW
        return m
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        return self

class _thin(_matrix):  ## the mother of all thin elements
    def __init__(self,particle=Particle.soll):
        self.particle = copy(particle)      ## keep a local copy of the particle instance (important!)
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
    particle=Particle.soll):
        super(QFth,self).__init__(particle=particle)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,particle=self.particle,label=self.label,viseo=+0.5)
        df = di
        kick = _matrix()    ## 6x6 unit matrix
        m = kick.matrix     ## thin lens quad matrix
        if(self.k0l == 0.):
#             m[1,0] = m[3,2] = 0.
            m[XPKOO,XKOO] = m[YPKOO,YKOO] = 0.
        else:
#             m[1,0] = -1./self.k0l
#             m[3,2] = -m[1,0]
            m[XPKOO,XKOO] = -1./self.k0l
            m[YPKOO,YKOO] = -m[XPKOO,XKOO]
        lens = df * (kick * di)     #matrix produkt df*kick*di
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = +0.5
    def shorten(self,l=0.):
        return self

class QDth(_thin):     ## thin D-quad
    def __init__(self,
    k0=0.,
    length=0.,
    label='QDT',
    particle=Particle.soll):
        super(QDth,self).__init__(particle=particle)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,particle=self.particle,label=self.label,viseo=-0.5)
        df = di
        kick = _matrix()    ## 6x6 unit matrix
        m = kick.matrix     ## thin lens quad matrix
        if(self.k0l == 0.):
#             m[1,0] = m[3,2] = 0.
            m[XPKOO,XKOO] = m[YPKOO,YKOO] = 0.
        else:
#             m[1,0] = 1./self.k0l
#             m[3,2] = -m[1,0]
            m[XPKOO,XKOO] = 1./self.k0l
            m[YPKOO,YKOO] = -m[XPKOO,XKOO]
        lens = (di * kick) * df
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return self

class RFC(_thin):      ## RF cavity as D*RFG*D
    def __init__(self,
                    U0=CONF['spalt_spannung'],
                    PhiSoll=radians(CONF['soll_phase']),
                    fRF=CONF['frequenz'],
                    label='RFC',
                    particle=Particle.soll,
                    gap=CONF['spalt_laenge'],
                    length=0.,
                    dWf=1.):
        super(RFC,self).__init__(particle=particle)
        self.u0     = U0*dWf
        self.phis   = PhiSoll
        self.freq   = fRF
        self.label  = label
        self.gap    = gap
        self.length = length
        self.dWf    = dWf
        self.di   = D(length=0.5*length,label='dcI',particle=self.particle)
        self.df   = D(length=0.5*length,label='dcF',particle=self.particle)
        self.kick = RFG(
                    U0=self.u0,
                    PhiSoll=self.phis,
                    fRF=self.freq,
                    label=self.label,
                    particle=self.particle,
                    gap=self.gap,
                    dWf=self.dWf)  ## Trace3D RF gap
#         DEBUG('RFC-kick deltaW >>',self.kick.deltaW)
        self.tr     = self.kick.tr
        tk_f = self.particle.tkin+self.kick.deltaW   #tkinetic after acc. gap
        self.df.adapt_for_energy(tk_f)               #update energy for downstream drift after gap
        lens = (self.df * self.kick) * self.di       #one for three
        self.matrix = lens.matrix
#         DEBUG('RFC matrix >>\n',self.matrix)
        self.triplet = (self.di,self.kick,self.df)
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
#         DEBUG('adapt_for_energy',tkin)
        self.__init__(
                    U0=self.u0,
                    PhiSoll=self.phis,
                    fRF=self.freq,
                    label=self.label,
                    particle=Particle(tkin=tkin),
                    gap=self.gap,
                    length=self.length,
                    dWf=self.dWf)
        return self

#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
class Test(_matrix):
    def __init__(self,a,b,c,d,e,f,label='test'):
        super(Test,self).__init__()
        self.matrix=NP.array([[a,b,0.,0.,0.,0.,0.,0.,0.,0.],
                              [c,d,0.,0.,0.,0.,0.,0.,0.,0.],
                              [0.,0.,a,b,0.,0.,0.,0.,0.,0.],
                              [0.,0.,d,e,0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,a,b,0.,0.,0.,0.],
                              [0.,0.,0.,0.,e,f,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                              ])
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
    print(a.string())
    b=Test(1,1,1,1,1,1,label='b')
    print(b.string())
    print((a*b).string())
    print((b*a).string())
    print('--------------- EOF test0 --------------------')
def test1():
    print('trivial test 1 ...')
    i1=_matrix()
    i2=i1*i1
    print(i1.string())
    print(i2.string())
    print('--------------- EOF test1 --------------------')
def test2():
    print('trivial test 2 ...')
    i1=_matrix()
    d1=D(10.,'D1')
    print(d1.string())
    print((d1*d1).string())
    print((d1*d1*d1).string())
    print((d1*i1*d1).string())
    d2=D(90,'D2')
    print(d2.string())
    print((d1*d2).string())
    d3=D(90.,label='')
    print((d2*d3).string())
    print('--------------- EOF test2 --------------------')
def test3():
    print('test product of _matrix class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.)
    print(qf.string())
    ## test product of _matrix class
    qd=QD(k0=k,length=1.)
    print(qd.string())
    print((qf*qd).string())
    print('--------------- EOF test3 --------------------')
def test4():
    print('test shortening of elements ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    ## elements
    d10=D(10.,'d10')
    print(d10.string())
    print((d10.shorten(1.e-2)).string())

    qf=QF(k0=k,length=1.)
    qf05=qf.shorten(0.2)
    print((qf05*qf05*qf05*qf05*qf05).string())
    print(qf.string())

    qd=QD(k0=k,length=1.)
    qd05=qd.shorten(0.2)
    print((qd05*qd05*qd05*qd05*qd05).string())
    print(qd.string())

    sd=SD(radius=10.,length=2.)
    sd05=sd.shorten(0.4)
    print((sd05*sd05*sd05*sd05*sd05).string())
    print(sd.string())
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
    print(mz.string())
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
            print(mi.string())
            m_end=m_end*mi
        print(m_end.string())
        print(m_anfang.string())
    print('--------------- EOF test6 --------------------')
def test7():
    print('======================================')
    print('test Rechteckmagnet...')
    rhob= wille()['beding_radius']
    lb=   wille()['dipole_length']
    mb=SD(radius=rhob,length=lb,label='B')
    mw=WD(mb,label='W')
    mr=mw*mb*mw
    print(mw.string())
    print(mb.string())
    print(mr.string())
    mr=RD(radius=rhob,length=lb,label='R')
    print(mr.string())
    print('--------------- EOF test7 --------------------')
def test8():
    print('test cavity...')
    objprnt(Particle.soll,'soll')
    cav=CAV()
    objprnt(cav,'CAV')
    objprnt(Particle.soll,'soll')
    rfg=RFG()
    objprnt(rfg,'RFG')
    objprnt(Particle.soll,'soll')
    print('--------------- EOF test8 --------------------')
def test9():
    print('\ntest: quad k-faktor and quad scaling')
    grad=CONF['quad_gradient']   # [T/m] gradient
    tk=CONF['injection_energy']    # [MeV]  kin. energy
    kq=k0(gradient=grad,tkin=tk)    # quad strength [1/m**2]
    len=0.4                         # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('\nproton {:.3f}[MeV] ~ beta {:.3f} in quadrupole:'.format(tk,Particle(tk).beta))
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
    Particle.soll = Particle(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqf.adapt_for_energy(tkf).string())
    print('========================')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    Particle.soll = Particle(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqd.adapt_for_energy(tkf).string())
    print('========================')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    Particle.soll = Particle(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(cavity.adapt_for_energy(tkf).string())
    print('--------------- EOF test9 --------------------')
def test10():
    print('\ntest: Particle class')
    dictprnt(CONF,text='setutils.CONF')

    # particle class
    print()
    print(Particle(0.).string())
    print(Particle(50.).string())
    print(Particle(200.).string())
    print(Particle(1.e6).string())
    print(Particle(1.e9).string())
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
    print(qf.string())
    print(qd.string())
    print('---------------- step through ---------------')
    for elm in qf.step_through(6):
        print(elm.string())
    for elm in qd.step_through(7):
        print(elm.string())
    print('------ RF cavity test & step through --------')
    objprnt(Particle.soll,'particle(i)')
    rf=RFC(length=length)
    print(rf.string())
    objprnt(Particle.soll,'particle(f)')
    objprnt(rf,'RFC cavity')
    for elm in rf.step_through():
        print(elm.string())
    print('--------------- EOF test11 --------------------')
def test12():
    print('\ntest12 adapt_for_energy change:')
#     d=D(length=99.); print('id >>',d); d.string()
#     d.adapt_for_energy(tkin=1000.); print('id >>',d); print(d.string())
#     qf=QF(k0=1.5,length=0.3); print('id >>',qf); print(qf.string())
#     qf.adapt_for_energy(tkin=200.); print('id >>',qf); print(qf.string())
#     qd=QD(k0=1.5,length=0.3); print('id >>',qd); print(qd.string())
#     qd.adapt_for_energy(tkin=200.); print('id >>',qd); print(qd.string())
    rfc=RFC(length=1.23); print('id >>',rfc); print(rfc.string())
    rfc.adapt_for_energy(tkin=200.); print('id >>',rfc); print(rfc.string())
if __name__ == '__main__':
    CONF['verbose']=3
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
    test12()
