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
import sys
from math import sqrt,sinh,cosh,sin,cos,fabs,tan,floor,modf,pi,radians,degrees,ceil
from copy import copy
import numpy as NP
# import warnings

from setutil import wille,PARAMS,FLAGS,dictprnt,objprnt,Proton,Electron,DEBUG,MarkerActions
from setutil import dBdxprot,scalek0prot,k0prot,I0,I1,arrprnt
from NamedObject import NamedObject
from ParamsObject import ParamsObject

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF
DEBUG_MAP    = DEBUG_OFF

## MDIM
MDIM=10        # dimension of matrices
# x        x'        y        y'        z       dp/p0     E        dE        s        l
XKOO = 0;XPKOO = 1;YKOO = 2;YPKOO = 3;ZKOO = 4;ZPKOO = 5;EKOO = 6;DEKOO = 7;SKOO = 8;LKOO = 9

NP.set_printoptions(linewidth=132,formatter={'float':'{:>8.5g}'.format})  #pretty printing

## the mother of all lattice elements (a.k.a. matrices)
class _Node(NamedObject,ParamsObject,object):
    """
    Base class for transfer matrices (owns its particle instance!)
    """
    # MDIMxMDIM matrices used here
    def __init__(self, particle=PARAMS['sollteilchen'], position=[0,0,0]):
        NamedObject.__init__(self)
        ParamsObject.__init__(self)
        self.matrix    = NP.eye(MDIM)     # MDIMxMDIM unit matrix
        self.particle  = copy(particle)   # keep a local copy of the particle instance (IMPORTANT!)
        self.position  = position         # [entrance,middle,exit]
        self.length    = 0.
        self['slice_min'] = 0.001            # default - minimal slice length
        self['viseo']     = 0                # default - invisible
    def __call__(self,n=MDIM,m=MDIM):
        return self.matrix[:n,:m]   ## return upper left nxm submatrix
    def string(self):
        n  = 42
        nx = 300
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
    def __mul__(self,other):
        product=NP.einsum('ij,jk',self.matrix,other.matrix)
        res=_Node()
        if (self.label == 'no label'):
            res.label=other.label
        else:
            res.label=self.label+'*'+other.label
        res.length=self.length+other.length
        res.matrix=product
        return res
    def reverse(self):
        raise RuntimeError('_Node:reverse not implemented!')
    def inverse(self):
        raise RuntimeError('_Node:inverse not implemented!')
        sys.exit(1)
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
        raise RuntimeError('FATAL: _Node.shorten(): virtual member function called!')
    def make_slices(self,anz=PARAMS['nbof_slices']):
        mr = None  ## ignore the very small rest
        slices = []

        if self.length == 0.:           ## zero length element (like WD or CAV)
            slices.append(self)

        else:
            step = self.length/anz         ## calc step size
            if step < self['slice_min']:
                step  = self['slice_min']

            (step_fraction_part,step_int_part) = modf(self.length/step)

            rest = step * step_fraction_part
            mx   = self.shorten(step)     # shorten element to step length
            if rest > 1.e-3:
                mr = self.shorten(rest)
            elif rest < 0.:
                raise RuntimeError('FATAL: negative resting step size when stepping through - STOP')
                sys.exit(1)


            for i in range(int(step_int_part)):
                slices.append(mx)
        if mr != None : slices += [mr]
        return slices
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
    def map(self,i_track):
        """
        Linear mapping of trjectory from (i) to (f)
        """
        f_track = self.matrix.dot(i_track)

        # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:
            f = f_track.copy()
            for i in range(len(f_track)-4):
                f[i] =f[i]*1.e3
            arrprnt(f,fmt='{:6.3g},',txt='mtx_map: ')

        return f_track
    def soll_map(self,i_track):
        """
        Linear mapping of trjectory from (i) to (f)
        """
        f_track = self.matrix.dot(i_track)
        return f_track
## unity matrix
class I(_Node):
    def __init__(self, label='I', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle, position=position)
        self.label    = label
## marker
class MRK(I):
    def __init__(self, label='MRK', particle=PARAMS['sollteilchen'], position=[0,0,0], actions=[]):
        super().__init__(particle=particle,position=position)
        self.label    = label
        self.actions  = actions
    def shorten(self,l=0):
        return self
    def do_actions(self):                # do actions attached to the marker
        for action in self.actions:
            MarkerActions[action]()
    def adjust_energy(self,tkin):
        self.__init__(label=self.label, particle=self.particle(tkin), position=self.position, actions=self.actions)
        return self
## Trace3D drift space
class D(I):
    """
    Trace3D drift space
    """
    def __init__(self, length=0., label='D', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle,position=position)
        self.label    = label
        self.length   = length
        m = self.matrix
        g = self.particle.gamma
        m[XKOO,XPKOO] = m[YKOO,YPKOO] = self.length
        m[ZKOO,ZPKOO] = self.length/(g*g)
        m[SKOO,LKOO]  = self.length     #delta-s
    def shorten(self,l=0.):
        return D(length=l,label=self.label,particle=self.particle)
    def adjust_energy(self,tkin):
        self.__init__(length=self.length,label=self.label,particle=self.particle(tkin),position=self.position)
        return self
## Trace3D focussing quad
class QF(D):
    """
    Trace3D focussing quad
    """
    def __init__(self, k0=0., length=0., label='QF', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle, position=position)
        self.label    = label
        self.length   = length
        self.k0       = k0         ## Quad strength [m**-2]
        self.matrix   = self._mx_()
        self['viseo'] = +0.5
    def shorten(self,l=0.):
        ret = QF(k0=self.k0, length=l, label=self.label, particle=self.particle, position=self.position)
        # DEBUG_MODULE('QF: ',self.__dict__)
        # DEBUG_MODULE('QF.shorten: ',ret.__dict__)
        return ret
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
            print('QF: neither QF nor QD! should never happen! - STOP')
            sys.exit(1)
        m[SKOO,LKOO]  = self.length     #delta-s
        return m
    def adjust_energy(self,tkin):
        ki = self.k0
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle, position=self.position)
        return self
## Trace3D defocusing quad
class QD(QF):
    """
    Trace3D defocussing quad
    """
    def __init__(self, k0=0., length=0., label='QD', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(k0=k0,length=length,label=label,particle=particle, position=position)
        self['viseo'] = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0, length=l, label=self.label, particle=self.particle, position=self.position)
## Trace3D sector bending dipole in x-plane
class SD(D):
    """
    Trace3d sector dipole in x-plane
    """
    def __init__(self, radius=0., length=0., label='SD', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle, position=position)
        self.label  = label
        self.length = length
        self.radius = radius
        self.matrix = self._mx_()
        self['viseo'] = 0.25
    def shorten(self,l=0.):
        return SD(radius=self.radius, length=l, label=self.label, particle=self.particle, position=self.position)
    def _mx_(self):
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
        m[SKOO,LKOO] = self.length     #delta-s
        return m
    def adjust_energy(self,tkin):
        ri = self.radius
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        rf = ri*cpf/cpi  # scale bending radius with new impulse
        self.__init__(radius=rf, length=self.length, label=self.label, particle=self.particle, position=self.position)
        return self
## Trace3D rectangular bending dipole in x-plane
class RD(SD):
    """
    Trace3D rectangular dipole x-plane
    """
    def __init__(self, radius=0., length=0., label='RD', particle=PARAMS['sollteilchen'],position=[0,0,0]):
        super().__init__(radius=radius, length=length, label=label, particle=particle, position=position)
        psi = 0.5*length/radius   # halber Kantenwinkel
        
        self.wd = _wedge(psi,radius,particle,position)  # wedge
        rd = self.wd * (self * self.wd)
        self.matrix = rd.matrix
    def make_slices(self,anz=PARAMS['nbof_slices']):
        # DEBUG_MODULEll('RD.make_slices: {} {:8.4f}'.format(self.label,self.length))
        sdshort = self.shorten(self.length/anz)
        slices = [self.wd]          # wedge @ entrance
        for i in range(anz):
            slices.append(sdshort)
        slices.append(self.wd)      # wedge @ exit
        # DEBUG_MODULE('slices',slices)
        return slices
## Trace3D wedge of rectangular bending dipole in x-plane
class _wedge(I):
    """
    Trace3d dipole wedge x-plane
    """
    def __init__(self, psi, radius, particle, position):
        super().__init__(particle=particle, position=position)
        self.label  = 'w'
        ckp = tan(psi)/radius
        m = self.matrix
        # MDIMxMDIM matrix
        # m[1,0]      = +ckp
        # m[3,2]      = -ckp
        m[XPKOO,XKOO] = +ckp
        m[YPKOO,YKOO] = -ckp
## zero length RF-gap nach Dr.Tiede & T.Wrangler (simple)
class GAP(D):
    """
    Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    NOTE: zu einfach: produziert keine long. dynamik wie Trace3D RFG!
    """
    def __init__(self,
                        U0         = PARAMS['spalt_spannung'],
                        PhiSoll    = radians(PARAMS['soll_phase']),
                        fRF        = PARAMS['frequenz'],
                        label      = 'GAP',
                        particle   = PARAMS['sollteilchen'],
                        gap        = PARAMS['spalt_laenge'],
                        position   = [0,0,0],
                        dWf        = FLAGS['dWf']):
        super().__init__(particle=particle, position=position)
        self.label  = label
        self.length = 0.
        self.u0     = U0                       # [MV] gap Voltage
        self.phis   = PhiSoll                  # [radians] soll phase
        self.freq   = fRF                      # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = PARAMS['lichtgeschwindigkeit']/self.freq# [m] RF wellenlaenge
        self.tr     = self._trtf_(self.particle.beta)         # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf      # T.Wrangler pp.221

        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        particle    = copy(self.particle)
        part_center = particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap center
        g           = part_center.gamma                       # gamma @ gap center
        # self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g)                  # transport matrix
        self['viseo']  = 0.25
    def _trtf_(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*self.freq*self.gap / (beta*PARAMS['lichtgeschwindigkeit'])
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
    def adjust_energy(self,tkin):
        self.__init__(
                    U0         = self.u0,
                    PhiSoll    = self.phis,
                    fRF        = self.freq,
                    label      = self.label,
                    particle   = self.particle(tkin),
                    gap        = self.gap,
                    position   = self.position,
                    dWf        = self.dWf)
        return self
## Basic RF-gap model from A.Shislo
class _rfb(D):
    """
    Base RF Gap Model from pyOrbit (A.Shislo)
    """
    def __init__(self,
                    parent,
                    particle   = PARAMS['sollteilchen']):
        super().__init__(particle=particle)
        self.label    = 'RFB'
        self.length   = 0.
        self.particle = particle
        self.parent   = parent
        self.u0       = parent.u0         # [MV] gap Voltage
        self.phis     = parent.phis       # [radians] soll phase
        self.tr       = parent.tr
        self.gap      = parent.gap
        self.deltaW   = parent.deltaW
        self.lamb     = parent.lamb
    def map(self,i_track,mapping):
        if mapping == 'simple':
            which_map = self.lin_map
        elif mapping == 'base':
            which_map = self.rfb_map
        else:
            raise RuntimeError('"map" enabled but wrong "mapping" for {} specified! - STOP'.format(self.parent.label))
            sys.exit(1)
        return which_map(i_track)
    def lin_map(self,i_track):
        """
        Mapping of track from position (i) to (f) in linear approx. (A.Shislo 4.1)
        """
        xi        = i_track[XKOO]       # [0]
        xpi       = i_track[XPKOO]      # [1]
        yi        = i_track[YKOO]       # [2]
        ypi       = i_track[YPKOO]      # [3]
        zi        = i_track[ZKOO]       # [4] z-z0
        zpi       = i_track[ZPKOO]      # [5] dp/p - (dp/p)0
        Ti        = i_track[EKOO]       # [6] summe aller delta-T
        si        = i_track[SKOO]       # [8] summe aller laengen

        qE0L       = self.u0
        m0c2       = self.particle.e0
        T          = self.tr
        lamb       = self.lamb
        phis       = self.phis
        qE0LT      = qE0L*T
        twopi      = 2.*pi

        particlesi = self.particle
        Wsi        = particlesi.tkin
        betasi     = particlesi.beta
        gammasi    = particlesi.gamma
        gbsi       = particlesi.gamma_beta

        DWs        = self.deltaW
        Wsf        = Wsi + DWs
        particlesf = copy(particlesi)(tkin=Wsf)
        betasf     = particlesf.beta
        gammasf    = particlesf.gamma
        gbsf       = particlesf.gamma_beta

        m11 = gbsf/gbsi
        m12 = 0.
        m21 = qE0LT*twopi/(lamb*betasi)*sin(phis)
        m22 = 1.
        condPdT = m0c2*betasi**2*gammasi
        DWi = condPdT*zpi  # dp/p --> dT

        # THE MAP (a 2x2 matrix which is always linear!) A.Shishlo (4.1.6-10)
        zf  = m11*zi + m12*DWi
        DWf = m21*zi + m22*DWi

        condTdP = 1./(m0c2*betasf**2*gammasf)
        zfp = DWf*condTdP   # dT --> dp/p

        xf   = xi     # x does not change
        yf   = yi     # y does not change
        Tf   = Ti+DWs
        sf   = si     # because self.length always 0
        xpf  = gbsi/gbsf*xpi - xi * (pi*qE0LT/(m0c2*lamb*gbsi*gbsi*gbsf)) * sin(phis) # A.Shishlo 4.1.11)
        ypf  = gbsi/gbsf*ypi - yi * (pi*qE0LT/(m0c2*lamb*gbsi*gbsi*gbsf)) * sin(phis)

        f_track = NP.array([xf,xpf,yf,ypf,zf,zfp,Tf,1.,sf,1.])

        # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:
            f = f_track.copy()
            for i in range(len(f_track)-4):
                f[i] =f[i]*1.e3
            arrprnt(f,fmt='{:6.3g},',txt='lin_map: ')

        return f_track
    def rfb_map(self,i_track):
        """
        Mapping of track from position (i) to (f) in Base-RF-Gap model approx. (A.Shislo 4.2)
        """
        # DEBUG_MODULE('i_track:\n',str(i_track))
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] summe aller delta-T
        s        = i_track[SKOO]       # [8] summe aller laengen

        ttf        = self.tr
        qE0LT      = self.u0*ttf
        m0c2       = self.particle.e0
        lamb       = self.lamb
        twopi      = 2.*pi

        soll      = self.particle
        PHIN      = self.phis                    # soll phi @ (i)
        WIN       = soll.tkin                    # soll energy @ (i)
        WOUT      = qE0LT*cos(PHIN) + WIN        # soll energy @ (f) (4.1.6) A.Shishlo
        betai     = soll.beta
        gammai    = soll.gamma
        gbi       = soll.gamma_beta

        win       = zp * m0c2*betai**2*gammai + WIN     # energy @ (i) dp/p --> dT
        pin       = -z * twopi/(betai*lamb)   + PHIN    # phase  @ (i)

        # particlei = copy(self.particle)(tkin=win)
        # betai     = particlei.beta
        # gbi       = particlei.gamma_beta

        r    = sqrt(x**2+y**2)                          # radial coordinate
        Kr   = (twopi*r)/(lamb*gbi)
        i0   = I0(Kr)                                   # bessel function I0
        i1   = I1(Kr)                                   # bessel function I1

        wout   = win + qE0LT*i0*cos(pin)           # energy @ (f)   (4.2.3) A.Shishlo

        # particlef = copy(self.particle)(tkin=wout)
        particlef = copy(self.particle)(tkin=WOUT)
        betaf     = particlef.beta
        gammaf    = particlef.gamma
        gbf       = particlef.gamma_beta

        dw     = wout- WOUT

        z      = betaf/betai*z                        # z @ (f)
        zpf    = 1./(m0c2*betaf**2*gammaf) * dw       # dW --> dp/p @ (f)

        T   = T + WOUT-WIN

        commonf = qE0LT/(m0c2*gbi*gbf)*i1             # common factor
        if r > 0.:
            xp  = gbi/gbf*xp - x/r*commonf*sin(pin)   # Formel 4.2.6 A.Shishlo
            yp  = gbi/gbf*yp - y/r*commonf*sin(pin)
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x,xp,y,yp,z,zpf,T,1.,s,1.])

        # for DEBUGGING
        if DEBUG_MAP == DEBUG_ON:
            f = f_track.copy()
            for i in range(len(f_track)-4):
                f[i] =f[i]*1.e3
            arrprnt(f,fmt='{:6.3g},',txt='rfb_map: ')

        return f_track
## Trace3D zero length RF-gap
class RFG(D):
    """
    Trace3D zero length Rf-gap
    """
    def __init__(self,
                    U0         = PARAMS['spalt_spannung'],
                    PhiSoll    = radians(PARAMS['soll_phase']),
                    fRF        = PARAMS['frequenz'],
                    label      = 'RFG',
                    particle   = PARAMS['sollteilchen'],
                    gap        = PARAMS['spalt_laenge'],
                    position   = [0,0,0],
                    mapping    = 'simple',
                    dWf        = FLAGS['dWf']):
        super().__init__(particle=particle,position=position)
        self.label   = label
        self.length  = 0.
        self.u0      = U0*dWf             # [MV] gap Voltage
        self.phis    = PhiSoll            # [radians] soll phase
        self.freq    = fRF                # [Hz]  RF frequenz
        self.gap     = gap
        self.dWf     = dWf
        self.mapping = mapping if FLAGS['map'] else 'T3D'
        self.lamb    = PARAMS['lichtgeschwindigkeit']/self.freq # [m] RF wellenlaenge
        self.tr      = self._trtf_(self.particle.beta)
        self.deltaW  = self.u0*self.tr*cos(self.phis)         # deltaW energy kick nach Trace3D
        self['viseo']= 0.25
        # DEBUG_MODULE('RFG: \n',self.particle.string())
        # DEBUG_MODULE('RFG: U0,phis,tr: {:8.4}, {:8.4}, {:8.4}'.format(self.u0,degrees(self.phis),self.tr))
        # DEBUG_MODULE('RFG: deltaW: {:8.6e}'.format(self.deltaW))
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        particle    = copy(self.particle)
        part_center = particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap cemter
        g           = part_center.gamma                       # gamma @ gap center
        particlei   = self.particle                           # particle @ (i)
        particlef   = particle(particlei.tkin+self.deltaW)    # particle @ (f)
        self.particlef = particlef
        # DEBUG_MODULE('RFG: beta i,c,f {:8.6f},{:8.6f},{:8.6f}'.format(particlei.beta,b,particlef.beta))
        # self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g,particlei,particlef)   # the LINEAR TRANSPORT matrix R

        # !!!!!  INSTANCIATE a MAP instead of using the R matrix
        if FLAGS['map']:
            self.rfb = _rfb(self,particle=self.particle)

    def _trtf_(self,beta):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = pi*self.freq*self.gap / PARAMS['lichtgeschwindigkeit']
        # DEBUG_MODULE('RFG: teta , beta>>',teta,beta)
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
    def adjust_energy(self,tkin):
        self.__init__(
                    U0         = self.u0,
                    PhiSoll    = self.phis,
                    fRF        = self.freq,
                    label      = self.label,
                    particle   = self.particle(tkin),
                    gap        = self.gap,
                    position   = self.position,
                    mapping    = self.mapping,
                    dWf        = self.dWf)
        return self
    def map(self,i_track):
        """
        Mapping of track from position (i) to (f)
        """
        if FLAGS['map']:
            # NOTE: mapping with _rfb-map
            f_track = self.rfb.map(i_track,self.mapping)
        else:
            # NOTE: linear mapping with T3D matrix
            f_track = super().map(i_track)
        return f_track
    def soll_map(self,i_track):
        f_track = super().soll_map(i_track)
        return f_track
## base of _thin Nodes
class _thin(_Node):
    """
    Base class for thin elements implemented as triplet D*Kick*D
    """
    def __init__(self,particle=PARAMS['sollteilchen'],position=[0,0,0]):
        super().__init__(particle=particle, position=position)
    def make_slices(self,anz=PARAMS['nbof_slices']):  # stepping routine through the triplet
        # DEBUG_MODULEll('_thin.make_slices: {} {:8.4f}'.format(self.label,self.length))
        anz1 = int(ceil(anz/2))
        di   = self.triplet[0]
        df   = self.triplet[2]
        kik  = self.triplet[1]
        d1   = di.shorten(di.length/anz1)
        d2   = df.shorten(df.length/anz1)
        slices = []
        for i in range(anz1):
            slices.append(d1)
        slices.append(kik)                  # the Kick
        for i in range(anz1):
            slices.append(d2)
        # DEBUG_MODULE('slices',slices)
        return slices
## thin F-quad
class QFth(_thin):
    """
    Thin F-Quad
    """
    def __init__(self, k0=0., length=0., label='QFT', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle,position=position)
        self.label     = label
        self.length    = length
        self.k0        = k0
        self['viseo']  = +0.5
        di = D(length=0.5*self.length,particle=self.particle)
        df = di
        kick = _kick(quad=self, particle=self.particle, position=position)    # MDIMxMDIM unit matrix
        lens = df * (kick * di)     #matrix produkt df*kick*di
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
    def adjust_energy(self,tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        ki = self.k0
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle, position=self.position)
        return self
##_kick
class _kick(I):
    def __init__(self, quad=None, particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(label='k',particle=particle)
        m = self.matrix                         # my thin lens quad matrix
        # m[1,0]      = -self.k0*L
        # m[3,2]      = -m[1,0]
        m[XPKOO,XKOO] = -quad.k0*quad.length
        m[YPKOO,YKOO] = -m[XPKOO,XKOO]
## thin D-quad
class QDth(QFth):
    """
    Thin D-Quad
    """
    def __init__(self, k0=0., length=0., label='QDT', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(k0=-k0,length=length,label=label,particle=particle, position=position)
        self['viseo']  = -0.5
## thin F-quad(x)
class QFthx(D):
    """
    Thin F-Quad   (express version of QFth)
    """
    def __init__(self, k0=0., length=0., label='QFT', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(particle=particle, position=position)
        self.label    = label
        self.length   = length
        self.k0       = k0
        self['viseo'] = +0.5
        L = self.length
        m = self.matrix                # thin lens quad matrix (by hand calculation)
        m[0,0]  = 1. - k0*(L**2)/2.
        m[0,1]  = L - k0*(L**3)/4.
        m[1,0]  = -k0*L
        m[1,1]  = m[0,0]
        m[2,2]  = 1. + k0*(L**2)/2.
        m[2,3]  = L + k0*(L**3)/4.
        m[3,2]  = +k0*L
        m[3,3]  = m[2,2]
    def adjust_energy(self,tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin)            # particle energy adjusted
        cpf = self.particle.gamma_beta
        ki  = self.k0
        kf  = ki*cpi/cpf               # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle, position=self.position)
        return self
    def make_slices(self,anz=PARAMS['nbof_slices']):
        slices = [self]
        return slices
## thin D-quad(x)
class QDthx(QFthx):
    """
    Thin D-Quad   (express version of QDth)
    """
    def __init__(self, k0=0., length=0., label='QDT', particle=PARAMS['sollteilchen'], position=[0,0,0]):
        super().__init__(k0=-k0,length=length,label=label,particle=particle, position=position)
        self['viseo'] = -0.5
## RF cavity als D*RFG*D
class RFC(_thin):
    """
    Rf cavity as product D*RFG*D (experimental!)
    """
    def __init__(self,
                    U0       = PARAMS['spalt_spannung'],
                    PhiSoll  = radians(PARAMS['soll_phase']),
                    fRF      = PARAMS['frequenz'],
                    label    = 'RFC',
                    particle = PARAMS['sollteilchen'],
                    gap      = PARAMS['spalt_laenge'],
                    length   = 0.,
                    position = [0,0,0],
                    dWf      = FLAGS['dWf']):
        super().__init__(particle=particle,position=position)
        if length == 0.: length = gap
        self.label  = label
        self.length = length
        self.u0     = U0*dWf
        self.phis   = PhiSoll
        self.freq   = fRF
        self.gap    = gap
        self.dWf    = dWf

        di   = D(length=0.5*length,particle=self.particle)
        df   = D(length=0.5*length,particle=self.particle)
        kick = RFG(
                    U0=self.u0,
                    PhiSoll=self.phis,
                    fRF=self.freq,
                    label=self.label,
                    particle=self.particle,
                    gap=self.gap,
                    dWf=self.dWf)  ## Trace3D RF gap
        self.tr = kick.tr
        tk_f = self.particle.tkin+kick.deltaW   #tkinetic after acc. gap
        df.adjust_energy(tk_f)                  #update energy for downstream drift after gap
        lens = df * kick * df         #one for three
        self.matrix = lens.matrix
        # DEBUG_MODULE('RFC matrix\n',self.matrix)
        self.triplet = (di,kick,df)
    def adjust_energy(self,tkin):
        self.__init__(
                    U0            = self.u0,
                    PhiSoll       = self.phis,
                    fRF           = self.freq,
                    label         = self.label,
                    particle      = self.particle(tkin),
                    gap           = self.gap,
                    length        = self.length,
                    position      = self.position,
                    dWf           = self.dWf)
        return self
## SixTrack drift map
class SIXD(D):
    """
    Drift with Sixtrack mapping (experimental!)
    """
    def __init__(self,length=0.,label="SIXD",particle=PARAMS['sollteilchen'],position=[0.,0.,0.]):
        super().__init__(length=length,particle=particle,position=position)
        self.label    = label
        self.length   = length
        self['viseo'] = 0.
        self.off_soll = copy(self.particle)
    def shorten(self,l=0.):
        return SIXD(length=l,label=self.label,particle=self.particle,position=self.position)
    def adjust_energy(self,tkin):
        self.__init__(length=self.length,label=self.label,particle=self.particle(tkin),position=self.position)
        return self
    def map(self,i_track):
        def fpsigma(psigma,soll):
            beta0    = soll.beta
            E0       = soll.e
            m0c2     = soll.e0
            res      = (1+beta0**2*psigma)**2-(m0c2/E0)**2
            res      = sqrt(res)/beta0-1.
            return res
        def einsplusfpsigma(psigma,soll):
            return 1.+fpsigma(psigma,soll)
        #conversion T3D ==> Ripken-Schmidt (sixtrack)
        def t3d2six(i_track):
            soll     = self.particle
            x        = i_track[XKOO]       # [0]
            xp       = i_track[XPKOO]      # [1]
            y        = i_track[YKOO]       # [2]
            yp       = i_track[YPKOO]      # [3]
            z        = i_track[ZKOO]       # [4] z
            dp2p     = i_track[ZPKOO]      # [5] dp/p
            T        = i_track[EKOO]       # [6] summe aller delta-T
            s        = i_track[SKOO]       # [8] summe aller laengen

            E0       = soll.e
            beta0    = soll.beta
            p0       = soll.p          # cp-soll [MeV]
            m0c2     = soll.e0
            p        = p0/(1.-dp2p)
            E        = sqrt(p**2+m0c2**2) #E aus dp2p und p0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            gb       = particle.gamma_beta
            beta     = particle.beta

            px       = gb*m0c2/E0*xp
            py       = gb*m0c2/E0*yp
            sigma    = z
            try:
                psigma = ((beta0/beta/(1.-dp2p))-1.)/beta0**2
            except:
                print('(dp2p,beta,beta0)',(dp2p,beta,beta0))
                print('in t3d2six(): bad psigma')
                sys.exit(1)
            f_track  = NP.array([x,px,y,py,sigma,psigma,T,1.,s,1.])
            return f_track
        # conversion Ripken-Schmidt (sixtrack) ==> T3D
        def six2t3d(i_track):
            soll     = self.particle
            x        = i_track[XKOO]
            px       = i_track[XPKOO]
            y        = i_track[YKOO]
            py       = i_track[YPKOO]
            sigma    = i_track[ZKOO]
            psigma   = i_track[ZPKOO]
            T        = i_track[EKOO]
            s        = i_track[SKOO]

            E0       = soll.e
            beta0    = soll.beta
            m0c2     = soll.e0
            eta      = beta0**2*psigma
            E        = (1.+eta)*E0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            beta     = particle.beta
            gb       = particle.gamma_beta

            xp       = px/(gb*m0c2/E0)
            yp       = py/(gb*m0c2/E0)
            z        = sigma
            dp2p     = 1.-beta0/beta/(1.+beta0**2*psigma)
            f_track  = NP.array([x,xp,y,yp,z,dp2p,T,1.,s,1.])
            return f_track
        # Ripken-Schmidt sixtrack map
        def rps_map(i_track,l):
            soll     = self.particle
            xi       = i_track[XKOO]
            pxi      = i_track[XPKOO]
            yi       = i_track[YKOO]
            pyi      = i_track[YPKOO]
            sigmai   = i_track[ZKOO]
            psigmai  = i_track[ZPKOO]
            T        = i_track[EKOO]
            s        = i_track[SKOO]

            E0       = soll.e
            beta0    = soll.beta
            m0c2     = soll.e0
            eta      = beta0**2*psigmai
            E        = (1.+eta)*E0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            beta     = particle.beta

            xf       = xi + pxi/einsplusfpsigma(psigmai,soll)*l
            pxf      = pxi
            yf       = yi + pyi/einsplusfpsigma(psigmai,soll)*l
            pyf      = pyi
            sigmaf   = sigmai + (1.-(beta0/beta)*(1.+0.5*(pxi**2+pyi**2)/einsplusfpsigma(psigmai,soll)**2))*l
            psigmaf  = psigmai
            f_track  = NP.array([xf,pxf,yf,pyf,sigmaf,psigmaf,T,1.,s,1.])
            return f_track
        ##body map
        f_track     = t3d2six(i_track)
        DEBUG_MODULE('t3d-->six\n',f_track)
        f_track     = rps_map(f_track,self.length)
        DEBUG_MODULE('SIXD.map\n',f_track)
        f_track     = six2t3d(f_track)
        DEBUG_MODULE('six-->t3d\n',f_track)

        f_track[SKOO] += self.length         # finally adjust total lattice length
        return f_track
## utilities
class Test(_Node):
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
        sys.exit(1)
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
    i1=_Node()
    i2=i1*i1
    print(i1.string())
    print(i2.string())
def test2():
    print('--------------------------------Test2---')
    print('trivial test 2 ...')
    i1=I()
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
    print('test product of _Node class ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient,beta,energy))
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    qf=QF(k0=k,length=1.1)
    print("QF-->",qf.string())
    # test product of _Node class
    qd=QD(k0=k,length=1.2)
    print("QD-->",qd.string())
    print("QF*QD-->",(qf*qd).string())
def test4():
    print('--------------------------------Test4---')
    def doit(elm,anz):
        elm_slices = elm.make_slices(anz=anz)
        print(''.join('{}\n'.format(el) for el in elm_slices))
        elmx = elm_slices[0]
        for slice in elm_slices[1:]:
            elmx = elmx*slice
        print(elmx.string())
        print(elm.string())
    print('test slicing of elements ...')
    gradient =1.
    beta     =0.5
    energy   =0.2
    k=k0test(gradient=gradient,energy=energy,beta=beta)
    # elements
    x=pi
    d       = D(length=x)
    qf      = QF(k0=k,length=x)
    qd      = QD(k0=k,length=x)
    sd      = SD(radius=10.,length=x)
    rd      = RD(radius=10.,length=x)
    gp      = GAP()
    rg      = RFG()
    qfth    = QFth(k0=k,length=x)
    qdth    = QDth(k0=k,length=x)
    qfthx   = QFthx(k0=k,length=x)
    qdthx   = QDthx(k0=k,length=x)
    rfc     = RFC(length=x)
    sixd    = SIXD(length=x)
    # slicing
    anz=5
    doit(d,anz)
    doit(sixd,anz)
    doit(qf,anz)
    doit(qd,anz)
    doit(sd,anz)
    doit(rd,anz)
    doit(gp,anz)
    doit(rg,anz)
    doit(qfth,anz)
    doit(qdth,anz)
    doit(qfthx,anz)
    doit(qdthx,anz)
    doit(rfc,anz)
def test5():
    print('--------------------------------Test5---')
    print("K.Wille's Beispiel auf pp.113 Formel (3.200)")
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
    mb  = RD(rhob,lb,'B')
    md  = D(ld)
    # test matrix multiplication
    mz = I()
    mz = mz *mqf
    mz = mz *md
    mz = mz *mb
    mz = mz *md
    mz = mz *mqd
    mz = mz *md
    mz = mz *mb
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
    mb=RD(rhob,lb,'B')
    md=D(ld)
    rfc=RFC(length=4*PARAMS['spalt_laenge'])

    steps = 13

    # test step_through elements ...
    list=[mqf,mqd,mb,md]
    list=[mqf,rfc]
    for m_anfang in list:
        m_end=I()
        slices = m_anfang.make_slices(anz = steps)
        print('======================================')
        for count,mi in enumerate(slices):
            print('step ',count+1,end='  ')
            print(mi.string())
            m_end=m_end*mi
        print(m_end,'\n'+m_end.string())
        print(m_anfang,'\n'+m_anfang.string())
def test7():
    print('--------------------------------Test7---')
    print('test Rechteckmagnet ...')
    rhob = wille()['bending_radius']
    lb   = wille()['dipole_length']
    mb   = SD(radius=rhob,length=lb,label='SD')
    mr   = RD(radius=rhob,length=lb,label='RD')
    print(mb.string())
    print(mr.string())
def test8():
    print('--------------------------------Test8---')
    print('test cavity ...')
    objprnt(PARAMS['sollteilchen'],'soll-particle')
    print('soll-particle\n'+PARAMS['sollteilchen'].string())
    cav = RFC()
    objprnt(cav,'RFC',filter='matrix')
    # objprnt(cav,'RFC')
    print('RFC.particle\n'+cav.particle.string())
    rfg = RFG()
    objprnt(rfg,'RFG',filter='matrix')
    # objprnt(rfg,'RFG')
    print('RFG.particle\n'+rfg.particle.string())
def test9():
    print('--------------------------------Test9---')
    print('test: quad k-faktor and quad scaling ...')
    grad = PARAMS['qd_gradient']         # [T/m] gradient
    tk   = PARAMS['injection_energy']    # [MeV]  kin. energy
    kq = k0prot(gradient=grad,tkin=tk) # quad strength [1/m**2]
    len=0.4                            # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('soll-particle\n'+PARAMS['sollteilchen'].string())
    print('kq [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))

    grad = dBdxprot(kq,tk) # quad gradient from k and tkinetic
    print('dB/dz[T/m]\t{:.3f} from dBxprot()'.format(grad))

    mqf = QF(k0=kq,length=len)
    mqd = QD(k0=kq,length=len)
    cavity = RFC(
        U0=PARAMS['spalt_spannung'],
        PhiSoll=radians(PARAMS['soll_phase']),
        fRF=PARAMS['frequenz'])
    print('======================== adjust_energy QF')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('(tki,kq) ({},{:.3f}) --> (tkf,k_scaled) ({},{:.3f})'.format(tki,kq,tkf,k_scaled))
        print(mqf.adjust_energy(tkf).string())
        print(mqf.particle.string())
    print('======================== adjust_energy QD')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('(tki,kq) ({},{:.3f}) --> (tkf,k_scaled) ({},{:.3f})'.format(tki,kq,tkf,k_scaled))
        print(mqd.adjust_energy(tkf).string())
        print(mqd.particle.string())
    print('======================== adjust_energy CAV')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('(tki,kq) ({},{:.3f}) --> (tkf,k_scaled) ({},{:.3f})'.format(tki,kq,tkf,k_scaled))
        print(cavity.adjust_energy(tkf).string())
        print(cavity.particle.string())
def test10():
    print('--------------------------------Test10---')
    print('Particle class test ...')
    dictprnt(PARAMS,text='setutil.PARAMS')
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
    print('thin lense tests ...')
    k0     = 1.
    length = 2.
    qf     = QFth(k0=k0,length=length)
    qd     = QDth(k0=k0,length=length)
    rf     = RFC(length=length)
    print(qf.string())
    print('soll-particle@QFT\n'+qf.particle.string())
    print(qd.string())
    print('soll-particle@QDT\n'+qf.particle.string())
    print(rf.string())
    print('soll-particle@RFC\n'+qf.particle.string())
    print('---------------- step through ---------------')
    for elm in qf.make_slices(anz=8):
        print(elm.string())
    print('---------------- step through ---------------')
    for elm in qd.make_slices(anz=7):
        print(elm.string())
    print('------ RF cavity test & step through --------')
    for elm in rf.make_slices():
        print(elm.string())
def test12():
    print('--------------------------------Test12---')
    print('test12 adjust_energy change ...')
    d = D(length=99.);           print('id >>',d);     print(d.string())
    d.adjust_energy(tkin=1000.); print('id >>',d);     print(d.string())
    qf = QF(k0=1.5,length=0.3);  print('id >>',qf);    print(qf.string())
    qf.adjust_energy(tkin=200.); print('id >>',qf);    print(qf.string())
    qd = QD(k0=1.5,length=0.3);  print('id >>',qd);    print(qd.string())
    qd.adjust_energy(tkin=200.); print('id >>',qd);    print(qd.string())
    rfc = RFC(length=1.23);      print('id >>',rfc);   print(rfc.string())
    rfc.adjust_energy(tkin=200.);print('id >>',rfc);   print(rfc.string())
def test13():
    print('--------------------------------Test13---')
    print('test SIXD node tracking ...')
    particle = Proton(tkin=100.)
    l    =  0.05    #[m]
    sixd = SIXD(length=l,particle=particle)
    xi   = yi  = 1.e-2
    xpi  = ypi = 1.e-2
    z    = 1.e-3
    dp2p = 1.e-2
    i_track = NP.array([xi,xpi,yi,ypi,z,dp2p,0.,1.,0.,1.])
    f_track = sixd.map(i_track)
    print(i_track)
    print(f_track)
def test14():
    print('--------------------------------Test14---')
    print('test MRO for QF,QD ...')
    print('type(QF.__mro__)= ',type(QF.__mro__))
    print(''.join('{}\n'.format(el) for el in QF.__mro__))
    print('type(QD.__mro__)= ',type(QD.__mro__))
    print(''.join('{}\n'.format(el) for el in QD.__mro__))
    
    qf0=QF(k0=1.0,length=1.98,label='qf0')
    qf1=QF(k0=1.0,length=1.88,label='qf1')
    qf1['viseo'] = 0.6
    qd0=QD(k0=1.0,length=1.78,label='qd0')
    qd0['viseo'] = -0.7
    print('qf0.label  = ',qf0.label)
    print('qf0["viseo"]= ',qf0['viseo'])
    print('qf1.label  = ',qf1.label)
    print('qf1["viseo"]= ',qf1['viseo'])
    print('qd0.label  = ',qd0.label)
    print('qd0["viseo"]= ',qd0['viseo'])
    print('\n',qf0.__dict__)
    print('\n',qd0.__dict__)
## main ----------
if __name__ == '__main__':
    FLAGS['verbose']=3
    # test0()
    # test1()
    # test2()
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
    test13()
    test14()
