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
import warnings

from setutil import wille,PARAMS,FLAGS,dictprnt,objprnt,Proton,Electron,DEBUG,MarkerActions
from setutil import dBdxprot,scalek0prot,k0prot,I0,I1,arrprnt

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
        self.slice_min = 0.001      ## minimal slice length
        self.viseo = 0.
    def __call__(self,n=MDIM,m=MDIM):
        return self.matrix[:n,:m]     ## return upper left nxm submatrix
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
        raise RuntimeError('FATAL: _matrix_.shorten(): virtual member function called!')
    def make_slices(self,anz=10):
        mb = I(label='',viseo=0)                   ## viseo point 
        mv = I(label='',viseo=self.viseo)          ## viseo point
        mr = I(label=self.label,viseo=self.viseo)  ## very small rest
        # slices = [mb,mv]
        slices = []

        if self.length == 0.:           ## zero length element (like WD or CAV)
            slices.append(self)
            # slices = slices + [mv,mb]
            # DEBUG_MODULE('_matrix_.make_slices',dict(label=self.label,length='zero lenth element'))

        else:
            step = self.length/anz         ## calc step size
            if step < self.slice_min:
                step  = self.slice_min

            (step_fraction_part,step_int_part) = modf(self.length/step)
            # DEBUG_MODULE('_matrix_.make_slices',dict(label=self.label,length=self.length,anz=anz,step=step,step_fracion=step_fraction_part,step_int=step_int_part))

            rest = step * step_fraction_part
            mx   = self.shorten(step)     # shorten element to step length
            if rest > 1.e-3:
                mr = self.shorten(rest)
            elif rest < 0.:
                raise RuntimeError('FATAL: negative resting step size when stepping through - STOP')
                sys.exit(1)
                

            for i in range(int(step_int_part)):
                slices.append(mx)
            # slices = slices + [mr,mv,mb]
        slices += [mr]
        # DEBUG_MODULE('_matrix_.make_slices',slices)
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
    def set_section(self,sec=''):
        """
        Setter for section tag (sections are not mandatory!)
        To distinguish different parts of the lattice, each element can be tagged by a section ID
        indicating the lattice part it belongs to.
        """
        self.sec = sec
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
## unity matrix (owns its particle instance!)
class I(_matrix_):     
    def __init__(self, label='I', viseo=0., particle=PARAMS['sollteilchen']):
        super().__init__()
        self.label=label
        self.viseo=viseo
        self.particle=copy(particle)  # keep a local copy of the particle instance (IMPORTANT!)
## marker
class MRK(I):        
    def __init__(self, label='MRK', particle=PARAMS['sollteilchen'], actions=[]):
        super().__init__(label=label, particle=particle)
        self.actions = actions
    def shorten(self,l=0):
        return self
    def do_actions(self):                   # do actions attached to the marker
        for action in self.actions:
            MarkerActions[action]()
    def adjust_energy(self,tkin):
        self.__init__(label=self.label, particle=self.particle(tkin), actions=self.actions)
        return self
## Trace3D drift space
class D(I):     
    """
    Trace3D drift space
    """
    def __init__(self, length=0., viseo=0., label='D', particle=PARAMS['sollteilchen']):
        super().__init__(label=label, viseo=viseo, particle=particle)
        self.length = length     ##  length [m]
        g = self.particle.gamma
        self.matrix[XKOO,XPKOO] = self.matrix[YKOO,YPKOO] = self.length
        self.matrix[ZKOO,ZPKOO] = self.length/(g*g)
        self.matrix[SKOO,LKOO]  = self.length     #delta-s
    def shorten(self,l=0.):    # returns a new instance!
        return D(length=l,label=self.label, particle=self.particle, viseo=self.viseo)
    def adjust_energy(self,tkin):
        self.__init__(length=self.length, viseo=self.viseo, label=self.label, particle=self.particle(tkin))
        return self
## Trace3D focussing quad
class QF(D):     
    """
    Trace3D focussing quad
    """
    def __init__(self, k0=0., length=0., label='QF', particle=PARAMS['sollteilchen']):
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
            print('QF: neither QF nor QD! should never happen! - STOP')
            sys.exit(1)
        return m
    def adjust_energy(self,tkin):
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
    def __init__(self, k0=0., length=0., label='QD', particle=PARAMS['sollteilchen']):
        super().__init__(k0=k0, length=length, label=label, particle=particle)
        self.viseo = -0.5
    def shorten(self,l=0.):
        return QD(k0=self.k0, length=l, label=self.label, particle=self.particle)
## Trace3D sector bending dipole in x-plane
class SD(D):         
    """
    Trace3d sector dipole in x-plane
    """
    def __init__(self, radius=0., length=0., label='SB', particle=PARAMS['sollteilchen']):
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
    def adjust_energy(self,tkin):
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
    def __init__(self, radius=0., length=0., label='RB', particle=PARAMS['sollteilchen']):
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
    def __init__(self, sector, label='WD', particle=PARAMS['sollteilchen']):
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
                        U0         = PARAMS['spalt_spannung'],
                        PhiSoll    = radians(PARAMS['soll_phase']),
                        fRF        = PARAMS['frequenz'],
                        label      = 'GAP',
                        particle   = PARAMS['sollteilchen'],
                        gap        = PARAMS['spalt_laenge'],
                        dWf        = FLAGS['dWf']):
        super().__init__(label=label, particle=particle)
        self.u0     = U0                       # [MV] gap Voltage
        self.phis   = PhiSoll                  # [radians] soll phase
        self.freq   = fRF                      # [Hz]  RF frequenz
        self.dWf    = dWf
        self.gap    = gap
        self.lamb   = PARAMS['lichtgeschwindigkeit']/self.freq  # [m] RF wellenlaenge
        self.tr     = self._trtf_(self.particle.beta)         # time-transition factor
        self.deltaW = self.u0*self.tr*cos(self.phis)*dWf      # T.Wrangler pp.221
        tk_center   = self.deltaW*0.5+self.particle.tkin      # energy in gap center
        particle    = copy(self.particle)
        part_center = particle(tk_center)                     # particle @ gap center
        b           = part_center.beta                        # beta @ gap center
        g           = part_center.gamma                       # gamma @ gap center
        # self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g)                  # transport matrix
        self.viseo  = 0.25
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
                    dWf        = self.dWf)
        return self
## Basic RF-gap model from A.Shislo
class RFB(D):
    """
    Base RF Gap Model from pyOrbit (A.Shislo)
    """
    def __init__(self,parent,
                    label      = 'RFB',
                    particle   = PARAMS['sollteilchen']):
        super().__init__(label=label, particle=particle)
        self.particle = particle
        self.label    = label
        self.parent   = parent
        self.viseo    = parent.viseo
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
        
        # DWs        = qE0LT*cos(phis) is same as self.deltaW
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

        f_track = NP.array([xf,xpf,yf,ypf,zf,zfp,Tf,1,sf,1])

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

        r    = sqrt(x**2+y**2)                                # radial coordinate
        Kr   = (twopi*r)/(lamb*gbi)
        i0   = I0(Kr)                                         # bessel function I0
        i1   = I1(Kr)                                         # bessel function I1

        wout   = win + qE0LT*i0*cos(pin)              # energy @ (f)   (4.2.3) A.Shishlo

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
                    mapping    = 'simple',
                    dWf        = FLAGS['dWf']):
        super().__init__(label=label, particle=particle)
        self.viseo   = 0.25
        self.u0      = U0*dWf                                 # [MV] gap Voltage
        self.phis    = PhiSoll                                # [radians] soll phase
        self.freq    = fRF                                    # [Hz]  RF frequenz
        self.label   = label
        self.gap     = gap
        self.dWf     = dWf
        self.mapping = mapping if FLAGS['map'] else 'T3D'
        self.lamb    = PARAMS['lichtgeschwindigkeit']/self.freq # [m] RF wellenlaenge
        self.tr      = self._trtf_(self.particle.beta)
        self.deltaW  = self.u0*self.tr*cos(self.phis)         # deltaW energy kick nach Trace3D
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
        # DEBUG_MODULE('RFG: beta i,c,f {:8.6f},{:8.6f},{:8.6f}'.format(particlei.beta,b,particlef.beta))
        # self.Ks     = 2.*pi/(self.lamb*g*b)                   # T.Wrangler pp.196
        self.matrix = self._mx_(self.tr,b,g,particlei,particlef)   # the LINEAR TRANSPORT matrix R
        self.particlei = particlei
        self.particlef = particlef
        # !!!!!  INSTANCIATE a MAP instead of using the R matrix
        if FLAGS['map']:
            self.rfb = RFB(self,label='RFB',particle=self.particle)

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
        # DEBUG_MODULE('RFG: adjust_energy',tkin)
        self.__init__(
                    U0         = self.u0,
                    PhiSoll    = self.phis,
                    fRF        = self.freq,
                    label      = self.label,
                    particle   = self.particle(tkin),
                    gap        = self.gap,
                    mapping    = self.mapping,
                    dWf        = self.dWf)
        return self
    def map(self,i_track):
        """
        Mapping of track from position (i) to (f)
        """
        if FLAGS['map']:
            # NOTE: mapping with RFB-map
            f_track = self.rfb.map(i_track,self.mapping)
        else:
            # NOTE: linear mapping with T3D matrix
            f_track = super().map(i_track)
        return f_track
    def soll_map(self,i_track):
        f_track = super().soll_map(i_track)
        return f_track
class _thin(_matrix_): 
    """
    Base class for thin elements implemented as triplet D*Kick*D
    """
    def __init__(self,particle=PARAMS['sollteilchen']):
        self.particle = copy(particle)     # keep a local copy of the particle instance (important!)
    def make_slices(self,anz=10):          # stepping routine through the triplet
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
        slices.append(kik)              ## the Kick
        for i in range(anz1):
            slices.append(d2)
        # DEBUG_MODULE('slices',slices)
        return slices
    def shorten(self,l=0.):
        warnings.warn("No need to shorten a thin element!")
        return self
    def set_section(self,sec=''):
        """
        Setter for section tag (sections are not mandatory!)
        To distinguish different parts of the lattice, each element can be tagged by a section ID
        indicating the lattice part it belongs to.
        """
        self.sec = sec
## thin F-quad
class QFthx(D):
    """
    Thin F-Quad   (express version of QFth)
    """
    def __init__(self, k0=0., length=0., label='QFT', particle=PARAMS['sollteilchen']):
        super().__init__(viseo=+0.5, length=length, label=label, particle=particle)
        self.k0     = k0
        self.length = length
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
        self.particle(tkin)          # particle energy adjusted
        cpf = self.particle(tkin).gamma_beta
        ki  = self.k0
        kf  = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle)
        return self
    def shorten(self,l=0.):
        return QFthx(k0=self.k0,length=l,label=self.label,particle=self.particle)
class QFth(_thin):   
    """
    Thin F-Quad
    """
    def __init__(self, k0=0., length=0., label='QFT', particle=PARAMS['sollteilchen']):
        super().__init__(particle=particle)
        self.k0     = k0
        self.length = length
        L = self.length
        self.label  = label
        di = D(length=0.5*length,particle=self.particle,label=self.label,viseo=+0.5)
        df = di
        kick = I(particle=self.particle)    # MDIMxMDIM unit matrix
        m = kick.matrix                     # thin lens quad matrix
            # m[1,0]      = -self.k0*L
            # m[3,2]      = -m[1,0]
        m[XPKOO,XKOO] = -self.k0*L
        m[YPKOO,YKOO] = -m[XPKOO,XKOO]
        lens = df * (kick * di)     #matrix produkt df*kick*di
        self.matrix = lens.matrix
        self.triplet = (di,kick,df)
        self.viseo = +0.5
    def adjust_energy(self,tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle(tkin).gamma_beta
        ki = self.k0
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle)
        return self
## thin D-quad
class QDthx(QFthx):
    """
    Thin D-Quad   (express version of QDth)
    """
    def __init__(self, k0=0., length=0., label='QDT', particle=PARAMS['sollteilchen']):
        super().__init__(k0 = -k0, length=length, label=label, particle=particle)
        self.viseo = -0.5
        self.k0 = k0   # hide parent's member
    def shorten(self,l=0.):
        return QDthx(k0=self.k0,length=l,label=self.label,particle=self.particle)
class QDth(QFth):
    """
    Thin D-Quad
    """
    def __init__(self, k0=0., length=0., label='QDT', particle=PARAMS['sollteilchen']):
        super().__init__(k0 = -k0, length=length, label=label, particle=particle)
        self.viseo = -0.5
        self.k0    = k0
    def adjust_energy(self,tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle(tkin).gamma_beta
        ki = self.k0
        kf = ki*cpi/cpf     # scale quad strength with new impulse
        self.__init__(k0=kf, length=self.length, label=self.label, particle=self.particle)
        return self
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
                    dWf      = FLAGS['dWf']):
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
        # DEBUG_MODULE('RFC-kick: deltaW: ',self.kick.deltaW)
        self.tr = self.kick.tr
        tk_f = self.particle.tkin+self.kick.deltaW   #tkinetic after acc. gap
        self.df.adjust_energy(tk_f)               #update energy for downstream drift after gap
        lens = self.df * self.kick * self.df         #one for three
        self.matrix = lens.matrix
        # DEBUG_MODULE('RFC matrix\n',self.matrix)
        self.triplet = (self.di,self.kick,self.df)
    def adjust_energy(self,tkin):
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
    rfc=RFC(length=4*PARAMS['spalt_laenge'])

    steps = 13

    # test step_through elements ...
    list=[mqf,mqd,mb,mw,md]
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
    objprnt(PARAMS['sollteilchen'],'soll')
    cav = CAV()
    objprnt(cav,'CAV')
    print('CAV.particle\n'+cav.particle.string())
    rfg = RFG()
    objprnt(rfg,'RFG')
    print('RFG.particle\n'+rfg.particle.string())
def test9():
    print('--------------------------------Test9---')
    print('test: quad k-faktor and quad scaling')
    grad = PARAMS['qd_gradient']         # [T/m] gradient
    tk   = PARAMS['injection_energy']    # [MeV]  kin. energy
    kq = k0prot(gradient=grad,tkin=tk) # quad strength [1/m**2]
    len=0.4                            # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('sollteilchen\n'+PARAMS['sollteilchen'].string())
    print('kq [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))

    grad = dBdxprot(kq,tk) # quad gradient from k and tkinetic
    print('dB/dz[T/m]\t{:.3f} from dBxprot()'.format(grad))

    mqf = QF(kq,len)
    mqd = QD(kq,len)
    cavity = CAV(
        U0=PARAMS['spalt_spannung'],
        PhiSoll=radians(PARAMS['soll_phase']),
        fRF=PARAMS['frequenz'])
    print('======================== adjust_energy QF')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqf.adjust_energy(tkf).string())
        print(mqf.particle.string())
    print('======================== adjust_energy QD')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(mqd.adjust_energy(tkf).string())
        print(mqd.particle.string())
    print('======================== adjust_energy CAV')
    tki=PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0.,950.]:
        tkf=tki+dt
        k_scaled = scalek0prot(kq,tki,tkf)
        print('k[{} MeV] {:.3f} --> k[{} MeV] {:.3f}'.format(tki,kq,tkf,k_scaled))
        print(cavity.adjust_energy(tkf).string())
        print(cavity.particle.string())
def test10():
    print('--------------------------------Test10---')
    print('Particle class test')
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
    print('test12 adjust_energy change:')
    d = D(length=99.);                            print('id >>',d);     print(d.string())
    d.adjust_energy(tkin=1000.);               print('id >>',d);     print(d.string())
    qf = QF(k0=1.5,length=0.3);                   print('id >>',qf);    print(qf.string())
    qf.adjust_energy(tkin=200.);               print('id >>',qf);    print(qf.string())
    qd = QD(k0=1.5,length=0.3);                   print('id >>',qd);    print(qd.string())
    qd.adjust_energy(tkin=200.);               print('id >>',qd);    print(qd.string())
    rfc = RFC(length=1.23);                       print('id >>',rfc);   print(rfc.string())
    rfc.adjust_energy(tkin=200.);              print('id >>',rfc);   print(rfc.string())
## main ----------
if __name__ == '__main__':
    FLAGS['verbose']=3
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
