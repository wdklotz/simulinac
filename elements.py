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
from math import sqrt,sinh,cosh,sin,cos,fabs,tan,floor,modf,pi,radians,degrees,ceil
from copy import copy
import numpy as NP
import warnings

from setutil import wille,CONF,dictprnt,objprnt,Proton,Electron,DEBUG
from setutil import dBdxprot,scalek0prot,k0prot

## MDIM
MDIM=10        # dimension of matrices
# x        x'        y        y'        z       dp/p0     E        dE        s        l     
XKOO = 0;XPKOO = 1;YKOO = 2;YPKOO = 3;ZKOO = 4;ZPKOO = 5;EKOO = 6;DEKOO = 7;SKOO = 8;LKOO = 9

NP.set_printoptions(linewidth=132,formatter={'float':'{:>8.5g}'.format})  #pretty printing

## the mother of all thick lattice elements (a.k.a. matrices)
class _matrix_(object):
    """
    Base class for thick matrices
    """
    # MDIMxMDIM matrices used here
    def __init__(self):
        self.matrix=NP.eye(MDIM)    ## MDIMxMDIM unit matrix
        self.label=''               ## default empty label
        self.length=0.              ## default zero length!
        self.slice_min = 0.005      ## minimal slice length
        self.viseo = 0.
    def string(self):
        n = 33
        nx = 200
        if len(self.label) > nx:
            label = self.label[:n]+'.....'+self.label[-n:]   ## when too long keep it short
        else:
            label = self.label
        try:
            s='{} [{}]\n'.format(label,self.sec)         ## sections are not mandatory
        except AttributeError:
            s='{}\n'.format(label)
        for i in range(MDIM):
            for j in range(MDIM):
                s+='{:8.4g} '.format(self.matrix[i,j])
            s+='\n'
        return s
    # def ppstring(self):
    #     import pprint
    #     pp = pprint.PrettyPrinter(indent=4)
    #     return pp.pformat(self)
    def __mul__(self,other):
        product=NP.einsum('ij,jk',self.matrix,other.matrix)
        res=_matrix_()
        if (self.label == ''):
            res.label=other.label
        else:
            res.label=self.label+'*'+other.label
        res.length=self.length+other.length
        res.matrix=product
        return res
    def reverse(self):
        raise RuntimeError('_matrix_:reverse not implemented!')
    def inverse(self):
        raise RuntimeError('_matrix_:inverse not implemented!')
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
        raise RuntimeError('FATAL: _matrix_.shorten(): virtual member function called!')
    def step_through(self,anz=10):
        """
        Generator function to step through an element:
        The central function to calculate s-dependent twiss functions (nontrivial!)
        Default is 10 steps/element.
        Minimal step size is self.slice_min.
        """
        mb = I(label='',viseo=0)                   ## viseo point 
        mv = I(label='',viseo=self.viseo)          ## viseo point
        slices = [mb,mv]

        if self.length == 0.:           ## zero length element (like WD or CAV)
            slices.append(self)
            slices = slices + [mv,mb]
            for slice in slices: yield slice
            return

        else:
            step = self.length/anz         ## calc step size
            if step < self.slice_min:
                step  = self.slice_min
            (step_fraction_part,step_int_part) = modf(self.length/step)

            rest = step * step_fraction_part
            mx   = self.shorten(step)     # shorten element to step length
            mr   = I(label=self.label,viseo=self.viseo)
            if rest > 1.e-3:
                mr = self.shorten(rest)
            elif rest < 0.:
                raise RuntimeError('FATAL: negative resting step size when stepping through')

            for i in range(int(step_int_part)):
                slices.append(mx)
            slices = slices + [mr,mv,mb]
            for slice in slices: yield slice
            return
    def beta_matrix(self):
        """
        The 6x6 matrix to track twiss functions through the lattice
        """
      # m11 =self.matrix[0,0];         m12 =self.matrix[0,1]
      # m21 =self.matrix[1,0];         m22 =self.matrix[1,1]
      # n11 =self.matrix[2,2];         n12 =self.matrix[2,3]
      # n21 =self.matrix[3,2];         n22 =self.matrix[3,3]
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
    def set_section(self,sec=''):
        """
        Setter for section tag (sections are not mandatory!)
        To distinguish different parts of the lattice, each element can be tagged by a section ID
        indicating the lattice part it belongs to.
        """
        self.sec = sec
## unity matrix (owns its particle instance!)
class I(_matrix_):     
    def __init__(self, label='I', viseo=0., particle=CONF['sollteilchen']):
        super().__init__()
        self.label=label
        self.viseo=viseo
        self.particle=copy(particle)  # keep a local copy of the particle instance (IMPORTANT!)
## marker
class MRK(I):        
    def __init__(self, label='MRK', particle=CONF['sollteilchen']):
        super().__init__(label=label, particle=particle)
    def shorten(self,l=0):
        return self
    def adapt_for_energy(self,tkin):
        self.__init__(label=self.label, particle=self.particle(tkin))
        return self
## Trace3D drift space
class D(I):     
    """
    Trace3D drift space
    """
    def __init__(self, length=0., viseo=0., label='D', particle=CONF['sollteilchen']):
        super().__init__(label=label, viseo=viseo, particle=particle)
        self.length = length     ##  length [m]
        g = self.particle.gamma
        self.matrix[XKOO,XPKOO] = self.matrix[YKOO,YPKOO] = self.length
        self.matrix[ZKOO,ZPKOO] = self.length/(g*g)
        self.matrix[SKOO,LKOO]  = self.length     #delta-s
    def shorten(self,l=0.):    # returns a new instance!
        return D(length=l,label=self.label, particle=self.particle, viseo=self.viseo)
    def adapt_for_energy(self,tkin):
        self.__init__(length=self.length, viseo=self.viseo, label=self.label, particle=self.particle(tkin))
        return self
## Trace3D focussing quad
class QF(D):     
    """
    Trace3D focussing quad
    """
    def __init__(self, k0=0., length=0., label='QF', particle=CONF['sollteilchen']):
        super().__init__(length=length, label=label, particle=particle)
        self.k0=k0         ## Quad strength [m**-2]
        self.matrix = self._mx_()
        self.viseo = +0.5
    def shorten(self,l=0.):
        return QF(k0=self.k0, length=l, label=self.label, particle=self.particle)
    def _mx_(self):
        m = self.matrix
        g = self.particle.gamma
        rzz12 = self.length/(g*g)
        kwurz = sqrt(self.k0)
        phi = self.length*kwurz
        # focusing
        cf  =cos(phi)
        sf  =sin(phi)/kwurz
        cfp =-kwurz*sin(phi)
        sfp =cf
        # defocusing
        cd = cosh(phi)
        sd = sinh(phi)/kwurz
        cdp = kwurz*sinh(phi)
        sdp = cd
        # MDIMxMDIM matrix
        if (isinstance(self,QF) and (isinstance(self,QD) == False)):
        #       0    0           0    1            1     0            1     1            2     2
            m[XKOO,XKOO]=cf; m[XKOO,XPKOO]=sf; m[XPKOO,XKOO]=cfp; m[XPKOO,XPKOO]=sfp; m[YKOO,YKOO]=cd
        #       2    3            3    2             3     3            4     5
            m[YKOO,YPKOO]=sd; m[YPKOO,YKOO]=cdp; m[YPKOO,YPKOO]=sdp; m[ZKOO,ZPKOO]=rzz12
        elif isinstance(self,QD):
            m[XKOO,XKOO]=cd; m[XKOO,XPKOO]=sd; m[XPKOO,XKOO]=cdp; m[XPKOO,XPKOO]=sdp; m[YKOO,YKOO]=cf; m[YKOO,YPKOO]=sf; m[YPKOO,YKOO]=cfp; m[YPKOO,YPKOO]=sfp; m[ZKOO,ZPKOO]=rzz12
        else:
            raise RuntimeError('QF._mx_: neither QF nor QD! should never happen!')
        return m
    def adapt_for_energy(self,tkin):
        ki = self.k0
        cpi = self.particle.gamma_beta
        cpf = self.particle(tkin).gamma_beta
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle(tkin))
        return self
## Trace3D defocusing
class QD(QF):       
    """
    Trace3D defocussing quad
    """
    def __init__(self, k0=0., length=0., label='QD', particle=CONF['sollteilchen']):
        super().__init__(k0=k0, length=length, label=label, particle=particle)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0, length=l, label=self.label, particle=self.particle)
## Trace3D sector bending dipole in x-plane
class SD(D):         
    """
    Trace3d sector dipole in x-plane
    """
    def __init__(self, radius=0., length=0., label='SB', particle=CONF['sollteilchen']):
        super().__init__(length=length, label=label, particle=particle)
        self.radius = radius
        self.matrix = self._mx_()
        self.viseo = 0.25
    def shorten(self,l=0.):
        return SD(radius=self.radius, length=l, label=self.label, particle=self.particle)
    def _mx_(self):  # nach Trace3D
        m = self.matrix
        rho = self.radius
        k = 1./rho
        phi = self.length/rho
        cx = cos(phi) ; sx = sin(phi)
        b = self.particle.beta
        # x-plane
        # m[0,0] = cx;          m[0,1] = sx/k;           m[0,5] = rho*(1.-cx)
        # m[1,0] = -sx*k;       m[1,1] = cx;             m[1,5] = sx
        m[XKOO,XKOO]  = cx;     m[XKOO,XPKOO]   = sx/k;  m[XKOO,ZPKOO]  = rho*(1.-cx)
        m[XPKOO,XKOO] = -sx*k;  m[XPKOO,XPKOO]  = cx;    m[XPKOO,ZPKOO] = sx
        # y-plane
        # m[2,3] = self.length
        m[YKOO,YPKOO] = self.length
        # z-plane
        # m[4,0] = -sx;       m[4,1] = -rho*(1.-cx);          m[4,5] = rho*sx-self.length*b*b
        m[ZKOO,XKOO] = -sx;   m[ZKOO,XPKOO] = -rho*(1.-cx);   m[ZKOO,ZPKOO] = rho*sx-self.length*b*b
        return m
    def adapt_for_energy(self,tkin):
        ri = self.radius
        cpi = self.particle.gamma_beta
        cpf = self.particle(tkin).gamma_beta
        rf = ri*cpf/cpi  # scale bending radius with new impulse
        self.__init__(radius=rf, length=self.length, viseo=self.viseo, label=self.label, particle=self.particle(tkin))
        return self
## Trace3D rectangular bending dipole in x-plane
class RD(SD):        
    """
    Trace3D rectangular dipole x-plane
    """
    def __init__(self, radius=0., length=0., label='RB', particle=CONF['sollteilchen']):
        super().__init__(radius=radius, length=length, label=label, particle=particle)
        wd = WD(self,label='',particle=particle)  # wedge myself...
        rd = wd * (self * wd)
        self.matrix= rd.matrix
    def shorten(self,l=0.):
        return RD(radius=self.radius, length=l, label=self.label, particle=self.particle)
## Trace3D wedge of rectangular bending dipole in x-plane
class WD(D):      
    """
    Trace3d dipole wedge x-plane
    """
    def __init__(self, sector, label='WD', particle=CONF['sollteilchen']):
        super().__init__(label=label, particle=particle)
        m = self.matrix
        self.parent = sector
        self.radius = sector.radius
        self.psi = sector.length/self.radius
        rinv = 1./self.radius
        psi = 0.5*self.psi  ## Kantenwinkel
        ckp = rinv*tan(psi)
        # MDIMxMDIM matrix
        # m[1,0]      = ckp
        # m[3,2]      = -ckp
        m[XPKOO,XKOO] = ckp
        m[YPKOO,YKOO] = -ckp
    def shorten(self,l=0.):
        wd = WD(self.parent, label=self.label, particle=self.particle)
        m = wd.matrix
        wd.psi = l/wd.radius
        rinv = 1./wd.radius
        psi = 0.5*wd.psi  ## Kantenwinkel
        ckp = rinv*tan(psi)
        # MDIMxMDIM matrix
        m[XPKOO,XKOO] = ckp
        m[YPKOO,YKOO] = -ckp
        return wd
## zero length RF-gap nach Dr.Tiede & T.Wrangler (simple)
class GAP(D):
    """
    Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    NOTE: zu einfach: produziert keine long. dynamik wie Trace3D RFG!
    """
    def __init__(self,
                        U0         = CONF['spalt_spannung'],
                        PhiSoll    = radians(CONF['soll_phase']),
                        fRF        = CONF['frequenz'],
                        label      = 'GAP',
                        particle   = CONF['sollteilchen'],
                        gap        = CONF['spalt_laenge'],
                        dWf        = 1.):
        super().__init__(label=label, particle=particle)
        self.u0     = U0                       # [MV] gap Voltage
        self.phis   = PhiSoll                  # [radians] soll phase
        self.freq   = fRF                      # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._trtf_(self.particle.beta)         # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf      # T.Wrangler pp.221
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        particle    = copy(self.particle)
        part_center = particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap center
        g           = part_center.gamma                       # gamma @ gap center
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g)                  # transport matrix
        self.viseo  = 0.25
    def _trtf_(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*self.freq*self.gap / (beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def _mx_(self,tr,b,g):               # cavity nach Dr.Tiede pp.33 
        m   = self.matrix
        e0  = self.particle.e0
        cyp = cxp = -pi*self.u0*tr*sin(self.phis)/(e0*self.lamb*g*g*g*b*b*b)  # T.Wrangler pp. 196
        # m[1,0]      = cxp
        # m[3,2]      = cyp
        # m[6,7]      = self.deltaW      # energy kick = acceleration
        m[XPKOO,XKOO] = cxp
        m[YPKOO,YKOO] = cyp
        m[EKOO,DEKOO] = self.deltaW      # energy kick = acceleration
        return m
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        self.__init__(
                    U0         = self.u0,
                    PhiSoll    = self.phis,
                    fRF        = self.freq,
                    label      = self.label,
                    particle   = self.particle(tkin),
                    gap        = self.gap,
                    dWf        = self.dWf)
        return self
## Trace3D zero length RF-gap
class RFG(D):       
    """
    Trace3D zero length Rf-gap
    """
    def __init__(self,
                    U0         = CONF['spalt_spannung'],
                    PhiSoll    = radians(CONF['soll_phase']),
                    fRF        = CONF['frequenz'],
                    label      = 'RFG',
                    particle   = CONF['sollteilchen'],
                    gap        = CONF['spalt_laenge'],
                    dWf        = 1.):
        super().__init__(label=label, particle=particle)
        self.viseo  = 0.25
        self.u0     = U0*dWf                                  # [MV] gap Voltage
        self.phis   = PhiSoll                                 # [radians] soll phase
        self.freq   = fRF                                     # [Hz]  RF frequenz
        self.label  = label
        self.gap    = gap
        self.dWf    = dWf
        self.lamb   = CONF['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._trtf_(self.particle.beta)
        self.deltaW = self.u0*self.tr*cos(self.phis)          # Trace3D
        # DEBUG('RFG: \n',self.particle.string())
        # DEBUG('RFG: U0,phis,tr: {:8.4}, {:8.4}, {:8.4}'.format(self.u0,degrees(self.phis),self.tr))
        # DEBUG('RFG: deltaW: {:8.6e}'.format(self.deltaW))
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        particle    = copy(self.particle)
        part_center = particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap cemter
        g           = part_center.gamma                       # gamma @ gap center
        particlei   = self.particle                           # particle @ entrance
        particlef   = particle(particlei.tkin+self.deltaW)    # particle @ exit
        # DEBUG('RFG: beta i,c,f {:8.6f},{:8.6f},{:8.6f}'.format(particlei.beta,b,particlef.beta))
        self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g,particlei,particlef)       # transport matrix
    def _trtf_(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = pi*self.freq*self.gap / CONF['lichtgeschwindigkeit']
        # DEBUG('RFG: teta , beta>>',teta,beta)
        teta = teta/beta
        ttf = sin(teta)/teta
        return ttf
    def _mx_(self,tr,b,g,particlei,particlef):   # RF gap nach Trace3D pp.17 (LA-UR-97-886)
        m = self.matrix
        e0 = self.particle.e0
        kz = 2.*pi*self.u0*tr*sin(self.phis)/(e0*b*b*self.lamb)
        ky = kx = -0.5*kz/(g*g)
        bgi = particlei.gamma_beta
        bgf = particlef.gamma_beta
        bgi2bgf = bgi/bgf
        # MDIMxMDIM
        # m[1,0]      = kx/bgf;    m[1,1]         = bgi2bgf
        # m[3,2]      = ky/bgf;    m[3,3]         = bgi2bgf
        # m[5,4]      = kz/bgf;    m[5,5]         = bgi2bgf
        m[XPKOO,XKOO] = kx/bgf;    m[XPKOO,XPKOO] = bgi2bgf
        m[YPKOO,YKOO] = ky/bgf;    m[YPKOO,YPKOO] = bgi2bgf
        m[ZPKOO,ZKOO] = kz/bgf;    m[ZPKOO,ZPKOO] = bgi2bgf
        #energy kick at acceleration
        # m[6,7]      = self.deltaW
        m[EKOO,DEKOO] = self.deltaW
        return m
    def shorten(self,l=0.):
        return self
    def adapt_for_energy(self,tkin):
        # DEBUG('RFG: adapt_for_energy',tkin)
        self.__init__(
                    U0         = self.u0,
                    PhiSoll    = self.phis,
                    fRF        = self.freq,
                    label      = self.label,
                    particle   = self.particle(tkin),
                    gap        = self.gap,
                    dWf        = self.dWf)
        return self
## the mother of all thin elements: keeps particle instance!
class _thin(_matrix_): 
    """
    Base class for thin elements
    """
    def __init__(self,particle=CONF['sollteilchen']):
        self.particle = copy(particle)      ## keep a local copy of the particle instance (important!)
    def step_through(self,anz=10):          ## stepping routine through the triplet (D,Kick,D)
        anz1 = int(ceil(anz/2))
        di   = self.triplet[0]
        df   = self.triplet[2]
        kik  = self.triplet[1]
        d1   = di.shorten(di.length/anz1)
        d2   = df.shorten(df.length/anz1)
        slices = []
        for i in range(anz1):
            slices.append(d1)
        slices.append(kik)              ## the Kick
        for i in range(anz1):
            slices.append(d2)

        for slice in slices:
            yield slice
    def shorten(self,l=0.):
        warnings.warn("No need to shorten a thin element!",RuntimeWarning)
        return self
    def set_section(self,sec=''):
        """
        Setter for section tag (sections are not mandatory!)
        To distinguish different parts of the lattice, each element can be tagged by a section ID
        indicating the lattice part it belongs to.
        """
        self.sec = sec

## thin F-quad
class QFth(_thin):   
    """
    Thin F-Quad
    """
    def __init__(self, k0=0., length=0., label='QFT', particle=CONF['sollteilchen']):
        super().__init__(particle=particle)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,particle=self.particle,label=self.label,viseo=+0.5)
        df = di
        kick = _matrix_()    ## MDIMxMDIM unit matrix
        m = kick.matrix      ## thin lens quad matrix
        if(self.k0l == 0.):
            # m[1,0]      = m[3,2] = 0.
            m[XPKOO,XKOO] = m[YPKOO,YKOO] = 0.
        else:
            # m[1,0]      = -1./self.k0l
            # m[3,2]      = -m[1,0]
            m[XPKOO,XKOO] = -1./self.k0l
            m[YPKOO,YKOO] = -m[XPKOO,XKOO]
        lens = df * (kick * di)     #matrix produkt df*kick*di
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = +0.5
    def adapt_for_energy(self,tkin):
        ki = self.k0
        cpi = self.particle.gamma_beta
        cpf = self.particle(tkin).gamma_beta
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle(tkin))
        return self
## thin D-quad
class QDth(_thin):  
    """
    Thin D-Quad
    """
    def __init__(self, k0=0., length=0., label='QDT', particle=CONF['sollteilchen']):
        super().__init__(particle=particle)
        self.k0     = k0
        self.length = length
        self.k0l    = k0*length
        self.label  = label
        di = D(length=0.5*length,particle=self.particle,label=self.label,viseo=-0.5)
        df = di
        kick = _matrix_()    ## MDIMxMDIM unit matrix
        m = kick.matrix      ## thin lens quad matrix
        if(self.k0l == 0.):
            # m[1,0]      = m[3,2]        = 0.
            m[XPKOO,XKOO] = m[YPKOO,YKOO] = 0.
        else:
            # m[1,0]      = 1./self.k0l
            # m[3,2]      = -m[1,0]
            m[XPKOO,XKOO] = 1./self.k0l
            m[YPKOO,YKOO] = -m[XPKOO,XKOO]
        lens = df * (kick * di)
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = -0.5
    def adapt_for_energy(self,tkin):
        ki = self.k0
        cpi = self.particle.gamma_beta
        cpf = self.particle(tkin).gamma_beta
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle(tkin))
        return self
## RF cavity als D*RFG*D
class RFC(_thin):    
    """
    Rf cavity as product D*RFG*D (experimental!)
    """
    def __init__(self,
                    U0       = CONF['spalt_spannung'],
                    PhiSoll  = radians(CONF['soll_phase']),
                    fRF      = CONF['frequenz'],
                    label    = 'RFC',
                    particle = CONF['sollteilchen'],
                    gap      = CONF['spalt_laenge'],
                    length   = 0.,
                    dWf      = 1.):
        super().__init__(particle=particle)
        if length == 0.: length = gap
        self.u0     = U0*dWf
        self.phis   = PhiSoll
        self.freq   = fRF
        self.label  = label
        self.gap    = gap
        self.length = length
        self.dWf    = dWf
        self.di   = D(length=0.5*length, label='dcI', particle=self.particle)
        self.df   = D(length=0.5*length, label='dcF', particle=self.particle)
        self.kick = RFG(
                    U0=self.u0,
                    PhiSoll=self.phis,
                    fRF=self.freq,
                    label=self.label,
                    particle=self.particle,
                    gap=self.gap,
                    dWf=self.dWf)  ## Trace3D RF gap
        # DEBUG('RFC-kick: deltaW: ',self.kick.deltaW)
        self.tr = self.kick.tr
        tk_f = self.particle.tkin+self.kick.deltaW   #tkinetic after acc. gap
        self.df.adapt_for_energy(tk_f)               #update energy for downstream drift after gap
        lens = self.df * (self.kick * self.di)       #one for three
        self.matrix = lens.matrix
        # DEBUG('RFC matrix\n',self.matrix)
        self.triplet = (self.di,self.kick,self.df)
    def adapt_for_energy(self,tkin):
        self.__init__(
                    U0            = self.u0,
                    PhiSoll       = self.phis,
                    fRF           = self.freq,
                    label         = self.label,
                    particle      = self.particle(tkin),
                    gap           = self.gap,
                    length        = self.length,
                    dWf           = self.dWf)
        return self
## utilities
class Test(_matrix_):
    def __init__(self,a,b,c,d,e,f,label='test'):
        super().__init__()
        self.matrix=NP.array([[ a, b,0.,0.,0.,0.,0.,0.,0.,0.],
                              [ c, d,0.,0.,0.,0.,0.,0.,0.,0.],
                              [0.,0., a, b,0.,0.,0.,0.,0.,0.],
                              [0.,0., d, e,0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0., a, b,0.,0.,0.,0.],
                              [0.,0.,0.,0., e, f,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,1.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,1.,0.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,1.,0.],
                              [0.,0.,0.,0.,0.,0.,0.,0.,0.,1.],
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
    print('--------------------------------Test0---')
    print('trivial test 0 ...')
    a=Test(1,2,3,4,5,6,label='a')
    print(a.string())
    b=Test(1,1,1,1,1,1,label='b')
    print(b.string())
    print((a*b).string())
    print((b*a).string())
def test1():
    print('--------------------------------Test1---')
    print('trivial test 1 ...')
    i1=_matrix_()
    i2=i1*i1
    print(i1.string())
    print(i2.string())
def test2():
    print('--------------------------------Test2---')
    print('trivial test 2 ...')
    i1=_matrix_()
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
def test3():
    print('--------------------------------Test3---')
    print('test product of _matrix_ class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.)
    print(qf.string())
    # test product of _matrix_ class
    qd=QD(k0=k,length=1.)
    print(qd.string())
    print((qf*qd).string())
def test4():
    print('--------------------------------Test4---')
    print('test shortening of elements ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    # elements
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
def test5():
    print('--------------------------------Test5---')
    print("K.Wille's Beispiel auf pp. 112-113")
    kqf  = wille()['k_quad_f']
    lqf  = wille()['length_quad_f']
    kqd  = wille()['k_quad_d']
    lqd  = wille()['length_quad_d']
    rhob = wille()['bending_radius']
    lb   = wille()['dipole_length']
    ld   = wille()['drift_length']
    # elements
    mqf = QF(kqf,lqf,'QF')
    mqd = QD(kqd,lqd,'QD')
    mb  = SD(rhob,lb,'B')
    mw  = WD(mb)
    md  = D(ld)
    # test matrix multiplication
    mz = I()
    mz = mz *mqf
    mz = mz *md
    mz = mz *mw
    mz = mz *mb
    mz = mz *mw
    mz = mz *md
    mz = mz *mqd
    mz = mz *md
    mz = mz *mw
    mz = mz *mb
    mz = mz *mw
    mz = mz *md
    mz = mz *mqf
    print(mz.string())
def test6():
    print('--------------------------------Test6---')
    print('test step_through elements ...')
    kqf=  wille()['k_quad_f']
    lqf=  wille()['length_quad_f']
    kqd=  wille()['k_quad_d']
    lqd=  wille()['length_quad_d']
    rhob= wille()['bending_radius']
    lb=   wille()['dipole_length']
    ld=   wille()['drift_length']

    # elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(mb)
    md=D(ld)
    rfc=RFC(length=4*CONF['spalt_laenge'])

    steps = 13

    # test step_through elements ...
    list=[mqf,mqd,mb,mw,md]
    list=[mqf,rfc]
    for m_anfang in list:
        m_end=I()
        print('======================================')
        for count,mi in enumerate(m_anfang.step_through(anz=steps)):
            print('step ',count+1,end='  ')
            print(mi.string())
            m_end=m_end*mi
        print(m_end,'\n'+m_end.string())
        print(m_anfang,'\n'+m_anfang.string())
def test7():
    print('--------------------------------Test7---')
    print('test Rechteckmagnet...')
    rhob = wille()['bending_radius']
    lb   = wille()['dipole_length']
    mb   = SD(radius=rhob,length=lb,label='B')
    mw   = WD(mb,label='W')
    mr=mw*mb*mw
    print(mw.string())
    print(mb.string())
    print(mr.string())
    mr = RD(radius=rhob,length=lb,label='R')
    print(mr.string())
def test8():
    print('--------------------------------Test8---')
    print('test cavity...')
    objprnt(CONF['sollteilchen'],'soll')
    cav = CAV()
    objprnt(cav,'CAV')
    print('CAV.particle\n'+cav.particle.string())
    rfg = RFG()
    objprnt(rfg,'RFG')
    print('RFG.particle\n'+rfg.particle.string())
def test9():
    print('--------------------------------Test9---')
    print('test: quad k-faktor and quad scaling')
    grad = CONF['qd_gradient']         # [T/m] gradient
    tk   = CONF['injection_energy']    # [MeV]  kin. energy
    kq = k0prot(gradient=grad,tkin=tk) # quad strength [1/m**2]
    len=0.4                            # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('sollteilchen\n'+CONF['sollteilchen'].string())
    print('kq [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))

    grad = dBdxprot(kq,tk) # quad gradient from k and tkinetic
    print('dB/dz[T/m]\t{:.3f} from dBxprot()'.format(grad))

    mqf = QF(kq,len)
    mqd = QD(kq,len)
    cavity = CAV(
        U0=CONF['spalt_spannung'],
        PhiSoll=radians(CONF['soll_phase']),
        fRF=CONF['frequenz'])
    print('======================== adapt_for_energy QF')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    CONF['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqf.adapt_for_energy(tkf).string())
        print(mqf.particle.string())
    print('======================== adapt_for_energy QD')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    CONF['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqd.adapt_for_energy(tkf).string())
        print(mqd.particle.string())
    print('======================== adapt_for_energy CAV')
    tki=CONF['injection_energy']    # [MeV]  kin. energy
    CONF['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(cavity.adapt_for_energy(tkf).string())
        print(cavity.particle.string())
def test10():
    print('--------------------------------Test10---')
    print('Particle class test')
    dictprnt(CONF,text='setutil.CONF')
    # particle class
    print()
    print( Proton(0.).string())
    print( Proton(50.).string())
    print( Proton(200.).string())
    print( Proton(1.e3).string())
    print()
    print( Electron(0.).string())
    print( Electron(50.).string())
    print( Electron(200.).string())
    print( Electron(1.e3).string())
def test11():
    print('--------------------------------Test11---')
    print('thin lense tests')
    k0     = 1.
    length = 2.
    qf     = QFth(k0=k0,length=length)
    qd     = QDth(k0=k0,length=length)
    rf     = RFC(length=length)
    print(qf.string())
    print('sollteilchen@QFthin\n'+qf.particle.string())
    print(qd.string())
    print('sollteilchen@QDthin\n'+qf.particle.string())
    print(rf.string())
    print('sollteilchen@RFC cavity\n'+qf.particle.string())
    print('---------------- step through ---------------')
    for elm in qf.step_through(6):
        print(elm.string())
    for elm in qd.step_through(7):
        print(elm.string())
    print('------ RF cavity test & step through --------')
    for elm in rf.step_through():
        print(elm.string())
def test12():
    print('--------------------------------Test12---')
    print('test12 adapt_for_energy change:')
    d = D(length=99.);                            print('id >>',d);     print(d.string())
    d.adapt_for_energy(tkin=1000.);               print('id >>',d);     print(d.string())
    qf = QF(k0=1.5,length=0.3);                   print('id >>',qf);    print(qf.string())
    qf.adapt_for_energy(tkin=200.);               print('id >>',qf);    print(qf.string())
    qd = QD(k0=1.5,length=0.3);                   print('id >>',qd);    print(qd.string())
    qd.adapt_for_energy(tkin=200.);               print('id >>',qd);    print(qd.string())
    rfc = RFC(length=1.23);                       print('id >>',rfc);   print(rfc.string())
    rfc.adapt_for_energy(tkin=200.);              print('id >>',rfc);   print(rfc.string())
## main ----------
if __name__ == '__main__':
    CONF['verbose']=3
    # test0()
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    test6()
    # test7()
    # test8()
    # test9()
    # test10()
    # test11()
    # test12()
