__version__='11.0.2.3'
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
# TODO adjust_energy() and shorten() return new objects. make sure object properties like links are correctly passed.
import sys
import unittest
from copy import copy
from separatrix import w2phi
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil import DEBUG_ON,DEBUG_OFF,Proton
from setutil import Ktp
from T3D_M   import T3D_G
import warnings
import math    as M
import numpy   as NP
import setutil as UTIL
import OXAL_M as OX

twopi = 2.*M.pi
# numpy pretty printing
NP.set_printoptions(linewidth = 132, formatter = {'float': '{:>8.5g}'.format})
def wrapRED(str): return UTIL.colors.RED+str+UTIL.colors.ENDC

""" ------- The mother of all lattice element objects (a.k.a. nodes)# ------ """
class Node(object):
    """ Base class for transfer matrices (linear map)
        ii)  is a dictionary (DictObject base class)
        ii)  each instance holds its copy of the refrence particle (self.particle)
    """
    def __init__(self):
        self.type         = self.__class__.__name__  # self's node type
        self.accelerating = False
        self.label        = None
        self.viseo        = None      # default - invisible
        self.particle     = None      # !!!IMPORTANT!!! local copy of the particle object
        self.matrix       = None      # MDIMxMDIM zero matrix used here
        self.position     = None      # [entrance, middle, exit]
        self.length       = None
        self.aperture     = None
        self.next         = None      # link to right Node
        self.prev         = None      # link to jeft Node
        self.twiss        = None      # twiss functions @ middle of Node
        self.sigxy        = None      # envelope function @ middle of Node
        self.ref_particle = None      # ref particle @ exit of Node
        self.ref_track    = None      # rerence track @ exit of Node
    def __call__(self, n = MDIM, m = MDIM):
        # return upper left n, m submatrix
        return self.matrix[:n, :m]
    def __mul__(self, other):
        """ define the (*) operator for Node objects """
        res = Node()
        if (self.label == '' or isinstance(self,DKD)):
            res.label = other.label
        else:
            res.label = self.label+'*'+other.label
        res.length = self.length + other.length
        """ matrix product """
        res.matrix = NP.dot(self.matrix, other.matrix)
        return res
    def map(self,i_track):
        """ standard mapping with T3D matrix """
        # ftrack = copy(itrack)    #TODO needed? think NO!
        f_track = NP.dot(self.matrix,i_track)
        return f_track
    def toString(self):
        ret = repr(self)
        for k,v in self.__dict__.items():
            ret+= '\n{}:{}'.format(k,v)
        return ret
    def adjust_energy(self, tkin):
        """ dummy adjust """
        return self
    def shorten(self, length):
        """ dummy shorten """
        return self
    def prmatrix(self):
        n  = 500
        nx = 1000
        if len(self.label) > nx:
            # make short when too long
            label = self.label[:n]+'.....'+self.label[-n:]
        else:
            label = self.label
        # try because sections are not mandatory
        try:
            s = '{} [{}]\n'.format(label, self.sec)
        except AttributeError:
            s = '{}\n'.format(label)
        for i in range(MDIM):
            for j in range(MDIM):
                s+='{:8.4g} '.format(self.matrix[i, j])
            s+='\n'
        return s
    def trace(self):
        return self.tracex()+self.tracey()
    def tracex(self):
        res = 0.
        for i in range(2):
            res += self.matrix[i, i]
        return res
    def tracey(self):
        res = 0.
        for i in range(2, 4):
            res += self.matrix[i, i]
        return res
    def make_slices(self, anz=1):
        slices = []
        if self.length == 0. or anz <=0:    anz = 1
        step = self.length/anz 
        mx   = self.shorten(step)
        for i in range(anz):
            slices.append(mx)
        return slices
    def beta_matrix(self):
        """ The 9x9 matrix to track twiss functions from the node's R-matrix """
        # Aliases
        m11  = self.matrix[XKOO, XKOO];   m12  = self.matrix[XKOO, XPKOO]
        m21  = self.matrix[XPKOO, XKOO];  m22  = self.matrix[XPKOO, XPKOO]
        
        n11  = self.matrix[YKOO, YKOO];   n12  = self.matrix[YKOO, YPKOO]
        n21  = self.matrix[YPKOO, YKOO];  n22  = self.matrix[YPKOO, YPKOO]

        o11  = self.matrix[ZKOO, ZKOO];   o12  = self.matrix[ZKOO, ZPKOO]
        o21  = self.matrix[ZPKOO, ZKOO];  o22  = self.matrix[ZPKOO, ZPKOO]

        m_beta  =  NP.array([
        [m11*m11,   -2.*m11*m12,       m12*m12,    0.,        0.,               0.,         0.,         0.,               0.],
        [-m11*m21,   m11*m22+m12*m21, -m22*m12,    0.,        0.,               0.,         0.,         0.,               0.],
        [m21*m21,   -2.*m22*m21,       m22*m22,    0.,        0.,               0.,         0.,         0.,               0.],
        [0.,        0.,                0.,         n11*n11,  -2.*n11*n12,       n12*n12,    0.,         0.,               0.],
        [0.,        0.,                0.,        -n11*n21,  n11*n22+n12*n21,  -n22*n12,    0.,         0.,               0.],
        [0.,        0.,                0.,         n21*n21,  -2.*n22*n21,       n22*n22,    0.,         0.,               0.],
        [0.,        0.,                0.,         0.,        0.,               0.,         o11*o11,   -2.*o11*o12,       o12*o12],
        [0.,        0.,                0.,         0.,        0.,               0.,        -o11*o21,    o11*o22+o12*o21, -o22*o12],
        [0.,        0.,                0.,         0.,        0.,               0.,         o21*o21,   -2.*o22*o21,       o22*o22]
        ])
        return m_beta
    def aper_check(self,new_tp,s,**kwargs):
        return False   
    @property
    def isAccelerating(self):
        return self.accelerating
class I(Node):
    """  Unity matrix: the unity Node """
    def __init__(self, label='I'):
        super().__init__()
        self.label    = label
        self.matrix   = NP.eye(MDIM,MDIM) # set the NODE's member variable
        self.length   = 0.
class MRK(I):
    """ 
    Marker node (a.k.a element): Each marker is parent of an agent that does the specific action.
    The action can be bypassed if the 'maction'-FLAG is False.
    """
    def __init__(self, label, active, viseo=1., particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.)):
        super().__init__(label)
        self.active     = active
        self.particle   = copy(particle)
        self.position   = position
        self.viseo      = viseo if self.active else 0.
    def no_action(self,*args):
        pass
class D(Node):
    """  Trace3D drift space  """
    def __init__(self, label, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), length=0.,aperture=None):
        super().__init__()
        self.label    = label
        self.particle = copy(particle)
        self.position = position
        self.length   = length
        self.aperture = aperture
        self.viseo    = 0.
        self.matrix   = NP.eye(MDIM,MDIM)
        m = self.matrix 
        g = self.particle.gamma
        m[XKOO, XPKOO] = m[YKOO, YPKOO] = self.length
        m[ZKOO, ZPKOO] = self.length/(g*g)
        m[SKOO, DSKOO]  = self.length # Ds the longitudinal length increase

    def adjust_energy(self, tkin):
        adjusted = D(self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortend =  D(self.label, particle=self.particle, position=(0.,0.,0.), length=length, aperture=self.aperture)
        return shortend
class DKD(D):
    """  Trace3D drift spacer for Drift-Kick-Drift sandwich (for use with DYNAC CAVNUM cavities) """
    def __init__(self, label, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), length=0.,aperture=None):
        super().__init__(label,particle=particle,position=position,length=length,aperture=aperture)

    def adjust_energy(self, tkin):
        adjusted = DKD(self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortend =  DKD(self.label, particle=self.particle, position=(0.,0.,0.), length=length, aperture=self.aperture)
        return shortend
class QF(Node):
    """ 
    Trace3D focussing quad  !!! NEUES API !!!!
    """
    def __init__(self, label, grad, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), length=0., aperture=None):
        super().__init__()
        self.label    = label
        self.grad     = abs(grad)                        # [T/m]
        self.particle = copy(particle)
        self.k02      = K(self.grad,self.particle)       # [1/m**2]
        self.position = position
        self.length   = length
        self.aperture = aperture
        self.viseo    = +0.5
        self.thin     = False
        self.matrix   = self._mx()
    def _mx(self):
        m = NP.eye(MDIM,MDIM)
        g = self.particle.gamma
        rzz12 = self.length/(g*g)
        k0w   = M.sqrt(self.k02)
        l = self.length
        k02 = self.k02
        phi = l*k0w
        f = 1./(k02*l)
        if f/l > 50.: self.thin = True
        DEBUG_OFF(F"{self.type}: k02={k02:.2f} \tf={f:.2f} >> l={l:.2f}? {self.thin}")
        if self.thin != True:
        # if True:
            """ thick quad """
            # focusing
            cf   = M.cos(phi)
            sf   = M.sin(phi)/k0w
            cfp  = -k0w*M.sin(phi)
            sfp  = cf 
            # defocusing
            cd  =  M.cosh(phi)
            sd  = M.sinh(phi)/k0w
            cdp = k0w*M.sinh(phi)
            sdp = cd
            if self.type == "QF":
                m[XKOO, XKOO]  = cf; m[XKOO, XPKOO] = sf; m[XPKOO, XKOO] = cfp; m[XPKOO, XPKOO] = sfp
                m[YKOO, YKOO]  = cd; m[YKOO, YPKOO] = sd; m[YPKOO, YKOO] = cdp; m[YPKOO, YPKOO] = sdp
            else:
                m[XKOO, XKOO]  = cd; m[XKOO, XPKOO] = sd; m[XPKOO, XKOO] = cdp; m[XPKOO, XPKOO] = sdp
                m[YKOO, YKOO]  = cf; m[YKOO, YPKOO] = sf; m[YPKOO, YKOO] = cfp; m[YPKOO, YPKOO] = sfp
        else:
            """ thin quad D-Kick-D"""
            if self.type == "QF":
                pass
            else:
                f = -f
            m[XKOO, XKOO]  = 1.-l/(2.*f); m[XKOO, XPKOO] = l-l**2/(6*f); m[XPKOO, XKOO] = -1./f; m[XPKOO, XPKOO] = m[XKOO, XKOO]
            f = -f
            m[YKOO, YKOO]  = 1.-l/(2.*f); m[YKOO, YPKOO] = l-l**2/(6*f); m[YPKOO, YKOO] = -1./f; m[YPKOO, YPKOO] = m[YKOO, YKOO]
        m[ZKOO, ZPKOO] = rzz12
        m[SKOO, DSKOO]  = self.length # length increase
        return m
    def isThin(self):
        return self.thin
    def adjust_energy(self, tkin):
        adjusted = QF(self.label, self.grad, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortened = QF(self.label, self.grad, particle=self.particle, position=self.position, length=length, aperture=self.aperture)
        return shortened
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()
        fifo_xy= kwargs['fifo_xy']
        sfifo_xy=kwargs['sfifo_xy']
        lost=False

        # transverse apertures
        if self.aperture != None and not (abs(new_point[Ktp.x]) < self.aperture or abs(new_point[Ktp.y]) < self.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost
class QD(QF):
    """ 
    Trace3D defocussing quad  !!! NEUES API !!!! 
    """
    def __init__(self, label, grad, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), length=0., aperture=None):
        super().__init__(label, grad, particle=particle, position=position, length=length, aperture=aperture)
        self.viseo = -0.5
    def adjust_energy(self, tkin):
        adjusted = QD(self.label, self.grad, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortened = QD(self.label, self.grad, particle=self.particle, position=self.position, length=length,  aperture=self.aperture)
        return shortened
class SD(Node):
    """ Trace3d horizontal sector magnet. n=0 pure dipole, alpha in [deg], rho in [m]."""
    def __init__(self, label, alpha, rho, n=0, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
        super().__init__()
        self.label    = label
        self.alpha    = alpha   # [deg]
        self.rho      = rho     # [m]
        self.n        = n
        self.particle = copy(particle)
        self.position = position
        self.length   = rho*M.radians(self.alpha)
        self.aperture = aperture
        self.viseo    = 0.25
        self.matrix   = self._mx()
    def _mx(self):
        m = NP.eye(MDIM,MDIM)
        beta  = self.particle.beta
        gamma = self.particle.gamma
        alpha = M.radians(self.alpha)      # [rad]
        rho   = self.rho                 # [m]
        Ds     = abs(rho*alpha)          # [m]
        h = alpha/Ds      # [1/m]
        kx = M.sqrt((1-self.n)*h**2)       # [1/m]
        ky = M.sqrt(self.n*h**2)
        self.length = Ds
        cx = M.cos(kx*Ds)
        sx = M.sin(kx*Ds)
        cy = M.cos(ky*Ds)
        sy = M.sin(ky*Ds)

        m[XKOO, XKOO]  = cx;      m[XKOO, XPKOO]    = sx/kx            # x,x'-plane
        m[XPKOO, XKOO] = -kx*sx;  m[XPKOO, XPKOO]   = m[XKOO, XKOO]        
        m[YKOO, YKOO]  = cy;      m[YKOO, YPKOO]    = sy/ky if self.n != 0. else self.length   # y,y'-plane
        m[YPKOO, YKOO] = -ky*sy;  m[YPKOO, YPKOO]   = m[YKOO, YKOO]
        # z,z'-plane
        m[XKOO,ZKOO]  = 0.;       m[XKOO,ZPKOO]   = h*(1.-cx)/kx**2    # Rxz
        m[XPKOO,ZKOO] = 0.;       m[XPKOO,ZPKOO]  = h*sx/kx
        m[ZKOO,XKOO]  = -h*sx/kx; m[ZKOO,XPKOO]   = -h*(1.-cx)/kx**2   # Rzx
        m[ZPKOO,XKOO] = 0.;       m[ZPKOO,XPKOO]  = 0.
        m[ZKOO,ZKOO]  = 1.;       m[ZKOO,ZPKOO]   = -1./rho**2*(kx*Ds*beta**2-sx)/kx**3+Ds*(1.-1./(rho**2*kx**2))/gamma**2   #Rzz
        m[ZPKOO,ZKOO] = 0.;       m[ZPKOO,ZPKOO]  = m[ZKOO,ZKOO]
        m[SKOO, DSKOO]  = self.length  # length increase
        return m
    def adjust_energy(self, tkin):
        adjusted = SD(self.label,self.alpha,self.rho,self.n,particle=self.particle(tkin),position=self.position,aperture=self.aperture)
        return adjusted
    def shorten(self, alpha):
        shortSD = SD(self.label,alpha,self.rho,self.n, particle=self.particle,position=self.position, aperture=self.aperture)
        return shortSD
    def make_slices(self, anz=2):
        shortSD = self.shorten(self.alpha/anz)
        slices = []
        for i in range(anz):
            slices.append(shortSD)
        return slices
class RD(SD):
    """ Trace3d horizontal rechteck magnet. n=0 pure dipole, alpha in [deg], rho in [m]."""
    def __init__(self, label, alpha, rho, wedge, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
        super().__init__(label, alpha, rho, n=0, particle=particle, position=position, aperture=aperture)
        self.wedge  = wedge
        self.matrix = NP.dot(self.wedge.matrix,NP.dot(self.matrix,self.wedge.matrix))
    def adjust_energy(self, tkin):
        adjusted = RD(self.label,self.alpha,self.rho,self.wedge,particle=self.particle(tkin),position=self.position,aperture=self.aperture)
        return adjusted
    def shorten(self, alpha):
        shortSD = RD(self.label,alpha,self.rho,self.wedge, particle=self.particle,position=self.position, aperture=self.aperture)
        return shortSD
    def make_slices(self, anz=2):
        if anz < 2: anz = 2
        slicewinkel   = self.alpha/anz
        sdi   = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
        sdi.matrix = NP.dot(self.wedge.matrix,sdi.matrix)
        slices  = [sdi]          # wedge @ entrance
        for i in range(1,anz-1):
            shortRD = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
            slices.append(shortRD)
        sdf   = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
        sdf.matrix = NP.dot(sdf.matrix,self.wedge.matrix)
        slices.append(sdf)
        return slices
class Wedge(Node):
    """  Trace3d kantenfokussierung .a.k.a wedge: Ryy simplified if t3d_wedge=False """
    def __init__(self, kwinkel, rho, t3d_wedge=True):
        super().__init__()
        self.kwinkel = kwinkel     # [deg ] kantenwinkel
        self.rho = rho
        self.t3d_wedge = t3d_wedge
        self.label = "W"
        self.length = 0.
        beta = M.radians(self.kwinkel)    
        g = 0.050    # fixed gap of 50 mm assumed
        K1=0.45
        K2=2.8
        psi  = K1*g/rho*((1+M.sin(beta)**2)/M.cos(beta))*(1-K1*K2*(g/rho)*M.tan(beta)) if self.t3d_wedge else 0.
        mxpx = M.tan(beta)/rho
        mypy = M.tan(beta-psi)/rho
        mx   = NP.eye(MDIM,MDIM)
        mx[XPKOO, XKOO] = mxpx
        mx[YPKOO, YKOO] = -mypy
        self.matrix = mx
class GAP(Node):
    """ Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    ... nicht sehr nuetzlich: produziert keine long. Dynamik wie Trace3D RFG!  """
    def __init__(self, label, EzAvg, phisoll, gap, freq, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=UTIL.FLAGS['dWf']):
        """ EzAvg [MV/m], phisoll [rad], gap [m], freq [Hz] """
        super().__init__()
        self.accelerating = True
        self.label    = label
        self.EzAvg    = EzAvg*dWf
        self.phisoll  = phisoll
        self.gap      = gap
        self.freq     = freq
        self.particle = copy(particle)
        self.position = position
        self.aperture = aperture
        self.dWf      = dWf
        self.viseo    = 0.25
        self.length   = 0.
        E0L           = self.EzAvg*self.gap            # Spaltspannung
        self.deltaW,self.matrix = self._mx(E0L)
    def _mx(self,E0L):
        """ cavity nach Dr.Tiede pp.33 """
        lamb   = UTIL.PARAMS['clight']/self.freq            # [m] wellenlaenge
        beta   = self.particle.beta                    # beta Einstein
        ttf    = self.ttf(beta,lamb,self.gap)          # time-transition factor
        bg     = self.particle.gamma_beta
        m0c2   = self.particle.e0
        deltaW = E0L*ttf*M.cos(self.phisoll)                   # delta-W T.Wrangler pp.221
        m      = NP.eye(MDIM,MDIM)
        cyp = cxp = -UTIL.pi*E0L*ttf*M.sin(self.phisoll)/(m0c2*lamb*bg**3)
        m[XPKOO, XKOO] = cxp
        m[YPKOO, YKOO] = cyp
        m[EKOO, DEKOO] = deltaW      # energy increase
        m[SKOO, DSKOO]  = self.length # length increase
        return deltaW,m
    def ttf(self, beta, lamb, gap):
        """ ttf-factor nach Panofsky (Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
        x = UTIL.pi*gap/(beta*lamb)
        ttf = NP.sinc(x/UTIL.pi)
        return ttf
    def adjust_energy(self, tkin):
        adjusted = GAP(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,particle=self.particle(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted
class RFG_OLD(Node):
    """ KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE """
    """  RF-gap of zero length with different kick gap-models """
    # def __init__(self, label, EzPeak, phisoll, gap, cavlen,freq, SFdata=0, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=UTIL.FLAGS['dWf'], mapping='t3d'):
    def __init__(self, label, **kwargs):
        super().__init__()
        self.label        = label
        self.length       = 0. # 0. because it's a kick
        self.viseo        = 0.25
        self.accelerating = True
        self.dWf          = UTIL.FLAGS['dWf']                 # dWf=1 with acceleration =0 else

        self.EzPeak    = kwargs.get('EzPeak',0)                 # [MV/m] peak gap field
        self.phisoll   = kwargs.get('phisoll',None)             # [radians] soll phase
        self.gap       = kwargs.get('gap',None)                 # [m] rf-gap
        self.cavlen    = kwargs.get('cavlen',None)              # [m] cavity length
        self.freq      = kwargs.get('freq',None)                # [Hz]  RF frequenz
        self.SFdata    = kwargs.get('SFdata',None)              # SuperFish data
        self.particle  = kwargs.get('particle',None)
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)
        self.mapping   = kwargs.get('mapping',None)             # map model

        self.EzPeak    = self.EzPeak*self.dWf         # [MV/m] peak gap field
        self.omega     = twopi*self.freq
        self.lamb      = UTIL.PARAMS['clight']/self.freq
        self.ttf       = None
        self.E0L       = None
        self.qE0LT     = None
        self.deltaW    = None
        self.particlef = None
        self.rf_gap    = None
    def dispatch_model_matrix(self):
        """ dispatching to different gap models """
        if self.mapping   == 't3d' :   #NOTE: t3d mapping is matrix multiplication
            self.matrix    = self.gap_object.T3D_matrix()
            # DEBUG_OFF('matrix',self.matrix)
            pass
        elif self.mapping == 'simple':
            self.matrix    = self.gap_object.simple_matrix()  #NOTE: simple mapping is matrix multiplication
            # DEBUG_OFF('matrix',self.matrix)
            pass
        elif self.mapping == 'oxal':
            (self.matrix,self.ttf,self.deltaW)  = self.gap_object.OXAL_matrix(self.particle.tkin) #NOTE: OXAL mapping is matrix multiplication
            DEBUG_ON(f'{self.gap_object},matrix',self.matrix)
            pass
        elif self.mapping == 'base':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            self.map = self.base_map_1
        elif self.mapping == 'ttf':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # TTF_G has its own mapping method
        elif self.mapping == 'dyn':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # DYN_G has its own mapping method
        # TODO other mappings not tested
        else:
            raise(UserWarning(f"INFO: RFG is a kick-model and does not work with {self.mapping} mapping! Use one of [t3d,simple,base,ttf,oxal]."))
            sys.exit()
    @property
    def gap_object(self):
        return self.rf_gap
    @gap_object.setter
    def gap_object(self,rf_gap):
        self.rf_gap = rf_gap
    def T3D_matrix(self):
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
            x = gap/(beta*lamb)
            res =NP.sinc(x)
            return res
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m              = NP.eye(MDIM,MDIM)
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.qE0LT     = self.E0L*self.ttf
        self.deltaW    = self.E0L*self.ttf*M.cos(self.phisoll)
        self.particlef = UTIL.Proton(self.particle.tkin+self.deltaW)
        Wavg    = self.particle.tkin+self.deltaW/2.   # average tkin
        pavg    = UTIL.Proton(Wavg)
        bavg    = pavg.beta
        gavg    = pavg.gamma
        m0c2    = pavg.e0
        kz      = twopi*self.E0L*self.ttf*M.sin(self.phisoll)/(m0c2*bavg*bavg*self.lamb)
        ky      = kx = -0.5*kz/(gavg*gavg)
        bgi     = self.particle.gamma_beta
        bgf     = self.particlef.gamma_beta
        bgi2bgf = bgi/bgf
        m       = NP.eye(MDIM,MDIM)
        m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
        m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
        m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
        m[EKOO, DEKOO] = self.deltaW
        m[SKOO, DSKOO]  = 0.
        return m
    def simple_matrix(self):
        """ Simplified Matrix Model. (A.Shislo 4.1) """
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
            x = gap/(beta*lamb)
            res =NP.sinc(x)
            return res
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.qE0LT     = self.E0L*self.ttf
        self.deltaW    = self.E0L*self.ttf*M.cos(self.phisoll)

        # xi   = i_track[XKOO]       # 0
        # xpi  = i_track[XPKOO]      # 1
        # yi   = i_track[YKOO]       # 2
        # ypi  = i_track[YPKOO]      # 3
        # zi   = i_track[ZKOO]       # 4 z
        # zpi  = i_track[ZPKOO]      # 5 Dp/p
        # T    = i_track[EKOO]       # 6 kinetic energy REF
        # dT   = i_track[DEKOO]      # 7 delta energy REF
        # S    = i_track[SKOO]       # 8 position REF
        # dS   = i_track[DSKOO]      # 9 delta position REF

        C         = UTIL.PARAMS['clight']
        particle  = self.particle
        m0c2      = particle.e0
        WIN       = particle.tkin
        betai     = particle.beta
        gammai    = particle.gamma
        gbi       = particle.gamma_beta

        deltaW    = self.deltaW
        lamb      = self.lamb
        qE0LT     = self.qE0LT
        phisoll   = self.phisoll

        WOUT       = WIN + deltaW
        particlef  = copy(particle)(tkin = WOUT)
        betaf      = particlef.beta
        gammaf     = particlef.gamma
        gbf        = particlef.gamma_beta

        condTdP = 1./(m0c2*betaf**2*gammaf)  # conversion DW --> Dp/p
        DDW = self.omega/C/betai*self.qE0LT*M.sin(self.phisoll) # A.Shishlo/J.Holmes (4.1.8)

        m = NP.eye(MDIM,MDIM)
        m[1,0] = m[3,2] = -(UTIL.pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf))*M.sin(phisoll)   # x',x = y',y   (4.1.11)
        # m[1,1] = m[3,3] = m[5,5] = gbi/gbf # x',x' = y',y'  m[5,5] like T3D
        m[1,1] = m[3,3] = gbi/gbf            # x',x' = y',y': m[5,5]=1 like in Shishlo's paper
        m[4,4] = betaf/betai                 # z,z
        m[5,4] = DDW*condTdP                 # z,z'
        m[6,7] = deltaW                      # T,dT
        return m
    def base_map_0(self, i_track):
        """alte map version bis 02.02.2022"""
        # def DEBUG_TRACK(inout,track):
        #     print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        """ Mapping (i) to (f) in Base RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy SOLL
        S        = i_track[SKOO]       # [8] position SOLL

        particle = self.particle
        m0c2     = particle.e0
        betai    = particle.beta
        gammai   = particle.gamma
        gbi      = particle.gamma_beta
        tki      = particle.tkin
        freq     = self.freq
        lamb     = self.lamb
        phisoll  = self.phisoll
        qE0LT    = self.qE0LT
        deltaW   = self.deltaW
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)
        max_r  = 0.05              # max radial excursion [m]
        r      = M.sqrt(x**2+y**2)  # radial coordinate
        if r > max_r:
            raise UTIL.OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = UTIL.I0(Kr)                               # bessel function I0
        i1     = UTIL.I1(Kr)                               # bessel function I1
        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)
        # SOLL
        WIN       = tki                               # energy (i)
        DELTAW    = deltaW                       # energy kick
        WOUT      = WIN + DELTAW                      # energy (f) (4.1.6) A.Shishlo/J.Holmes
        # PARTICLE
        converter = UTIL.WConverter(WIN,freq)
        # phin      = -z * twopi/(betai*lamb) + phis     # phase (i)  alte methode
        phin      = converter.zToDphi(z) + phisoll          # phase (i)
        deltaW    = qE0LT*i0*M.cos(phin)                   # energy kick
        # win     = (zp * (gammai+1.)/gammai +1.) * WIN  # energy (i) dp/p --> dT alte methode
        win       =  converter.Dp2pToDW(zp) + WIN         # energy (i) dp/p --> dT
        wout      = win + deltaW                         # energy (f)   (4.2.3) A.Shishlo/J.Holmes
        dw        = wout - WOUT                          # d(deltaW)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particlef = UTIL.Proton(WOUT)       # !!!IMPORTANT!!! SOLL particle (f)
        betaf     = particlef.beta
        gammaf    = particlef.gamma
        gbf       = particlef.gamma_beta

        # converter = WConverter(WOUT,frq)
        z         = betaf/betai*z                     # z (f) (4.2.5) A.Shishlo/J.Holmes
        # zpf     = gammaf/(gammaf+1.) * dw/WOUT      # dW --> dp/p (f)  alte methode
        zpf       = converter.DWToDp2p(dw)            # dW --> dp/p (f)
        # if 0: print('z ',z,'zpf ',zpf)

        commonf = qE0LT/(m0c2*gbi*gbf)*i1             # common factor
        if r > 0.:
            xp  = gbi/gbf*xp - x/r*commonf*M.sin(phin)  # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbf*yp - y/r*commonf*M.sin(phin)  # should be phi-middle
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x, xp, y, yp, z, zpf, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        if 0:
            UTIL.arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
            UTIL.arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        # the parent reads these attributes below
        self.particlef = particlef
        return f_track
    def base_map_1(self, i_track):
        """Neue map Version ab 03.02.2022 ist ein Remake um Korrecktheit der Rechnung zu testen. 
           Produziert dasselbe Verhalten wie base_map_0 """
        # def DEBUG_TRACK(inout,track):
        #     print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        """ Mapping (i) to (O) in Base RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy ref Teilchen
        S        = i_track[SKOO]       # [8] position gap

        particleRi = self.particle   # ref Teilchen (I)
        m0c2       = particleRi.e0
        betai      = particleRi.beta
        gammai     = particleRi.gamma
        gbi        = particleRi.gamma_beta
        wRi        = particleRi.tkin
        freq       = self.freq
        lamb       = self.lamb
        phisoll    = self.phisoll
        deg_phisoll= M.degrees(phisoll)
        qE0LT      = self.qE0LT
        deltaW     = self.deltaW
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)

        max_r  = 0.05              # max radial excursion [m]
        r      = M.sqrt(x**2+y**2)   # radial coordinate
        if r > max_r:
            raise UTIL.OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = UTIL.I0(Kr)            # bessel function I0
        i1     = UTIL.I1(Kr)            # bessel function I1

        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)

        # ref Teilchen
        wRo = wRi + deltaW                           # ref Teilchen energy (O)
 
        # Teilchen
        converter   = UTIL.WConverter(wRi,freq)
        deg_converter = M.degrees(converter.zToDphi(z)) 
        phiin       = converter.zToDphi(z) + phisoll 
        deg_phiin   = M.degrees(phiin)        # Teilchen phase (I)
        wo_wi       = qE0LT*i0*M.cos(phiin)                 # energy kick (Shislo 4.2.3)
        wi          =  converter.Dp2pToDW(zp) + wRi        # Teilchen energy (I) dp/p --> dT
        wo          = wi + wo_wi                          # Teilchen energy (O)   
        dw          = wo - wRo                            # Differenz der energy kicks von Teilchen und ref Teilchen (entspricht delta**2)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particleRo = UTIL.Proton(wRo)
        betao      = particleRo.beta
        gammao     = particleRo.gamma
        gbo        = particleRo.gamma_beta

        zo         = betao/betai*z                     # z (O) (4.2.5) A.Shishlo/J.Holmes
        zpo        = converter.DWToDp2p(dw)            # dW --> dp/p (O)

        # if 0: print('z ',z,'zpf ',zpf)

        factor = qE0LT/(m0c2*gbi*gbo)*i1               # common factor
        if r > 0.:
            xp  = gbi/gbo*xp - x/r*factor*M.sin(phiin)   # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbo*yp - y/r*factor*M.sin(phiin)
        elif r == 0.:
            xp  = gbi/gbo*xp
            yp  = gbi/gbo*yp

        f_track = NP.array([x, xp, y, yp, zo, zpo, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        # if 0:
        #     arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        #     arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')

        # """ the parent reads these attributes below """
        self.particlef = particleRo
        return f_track
    def adjust_energy(self, tkin):
        if self.gap_object == self:
            adjusted = RFG(self.label,EzPeak=self.EzPeak,phisoll=self.phisoll,gap=self.gap,cavlen=self.cavlen,freq=self.freq,SFdata=self.SFdata,particle=UTIL.Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf,mapping=self.mapping)
            adjusted.gap_object = adjusted
            adjusted.dispatch_model_matrix()
        elif self.mapping == 'oxal':
            self.particle = UTIL.Proton(tkin)
            self.gap_object.OXAL_matrix(tkin)
            adjusted = self
            # adjusted.dispatch_model_matrix()
        return adjusted
    def waccept(self):
        """ 
        Calculate longitudinal acceptance, i.e. phase space ellipse parameters: T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0*Dphi0 = ellipse_area/pi
        """
        rf_gap    = self.gap_object      # this RF gap to use: can be self or others like OXAL_G or TTF_G

        Ez0       = rf_gap.EzPeak
        ttf       = rf_gap.ttf
        phisoll   = rf_gap.phisoll         # [rad]
        lamb      = rf_gap.lamb            # [m]
        freq      = rf_gap.freq            # [Hz]
        particle  = rf_gap.particle

        E0T       = Ez0*ttf              # [MV/m]
        m0c2      = particle.e0          # [MeV]
        gb        = particle.gamma_beta
        beta      = particle.beta
        gamma     = particle.gamma
        tkin      = particle.tkin
        # DEBUG_OFF("waccept",dict(E0T=E0T,phisoll=degrees(phisoll),lamb=lamb,freq=freq,m0c2=m0c2,gb=gb,beta=beta,gamma=gamma,tkin=tkin))

        # converter for this node
        conv = UTIL.WConverter(tkin,freq)

        try:
            # LARGE amplitude oscillations (T.Wangler pp. 175 6.28). w = Dgamma = DW/m0c2 normalized energy spread """
            # DEBUG_OFF(f'w2phi {(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll)}')                                                                                                                                                              
            w0large = M.sqrt(w2phi(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll))
            # DEBUG_OFF(f'w0large {w0large}')                                                                                                                                                              
        except ValueError as ex:
            exception = ex
            w0large = -1
        try:
            # SMALL amplitude oscillations separatrix (T.Wangler pp.185) """
            w0small = M.sqrt(2.*E0T*gb**3*lamb*phisoll**2*M.sin(-phisoll)/(M.pi*m0c2))
            # DEBUG_OFF(f'w0small {w0small}')                                                                                                                                                              
        except ValueError as ex:
            exception = ex
            w0small = -1
        if w0large != -1: 
            wmax = w0large
        elif w0large == -1 and w0small != -1:
            wmax = w0small
        else:
            DEBUG_ON(f'{exception} reason: ttf={rf_gap.ttf}, E0T={E0T}')
            sys.exit(1)

        # this node Dp/p max on separatrix
        Dp2pmax = conv.wToDp2p(wmax) 

        #  convert T.Wangler units {Dphi,w} to {z,dp/p} units with 1st cavity parameters
        betaw_i,alfaw_i,gammaw,emitw_i = UTIL.PARAMS['twiss_w_i']()
        Dphi0_i = UTIL.PARAMS['Dphi0_i']
        w0_i = (gamma-1.)*UTIL.PARAMS['DT2T_i']
        z0_i,Dp2p0_i,emitz_i,betaz_i = conv.wtoz((Dphi0_i,w0_i,emitw_i,betaw_i))
        alfaz_i = 0.

        # omega sync for this node
        omgl0zuomg = M.sqrt(E0T*lamb*M.sin(-phisoll)/(2*M.pi*m0c2*gamma**3*beta))
        omgl_0     = omgl0zuomg*twopi*freq   # [Hz]

        # longitudinal acceptance check (always done)     #TODO  is this needed?
        # if wmax <= w0_i:
        #     si,sm,sf = self.position
        #     warnings.showwarning(
        #         colors.RED+'out of energy acceptance @ s={:.1f} [m]'.format(si)+colors.ENDC,
        #         UserWarning,'elements.py',
        #         'waccept')

        # phase acceptance (REMARK: phase limits are not dependent on Dp/p aka w)
        phi_2=2.*phisoll
        phi_1=-phisoll

        res =  dict (
                emitw_i         = emitw_i,      # 1st cavity emittance {Dphi,w} units [rad,1]
                z0_i            = z0_i,         # 1st cavity ellipse z-axe crossing (1/2 axis) [m]
                Dp2p0_i         = Dp2p0_i,      # 1st cavity ellipse dp/p-axe crossing (1/2 axis)
                twiss_z_i       = UTIL.Twiss(betaz_i, alfaz_i, emitz_i), # 1st cavity twis parameters
                DWmax           = wmax*m0c2,    # this node max delta-W on separatrix [MeV]
                Dp2pmax         = Dp2pmax,      # this node Dp/p max on separatrix [1]
                phaseacc        = (conv,phi_2,phisoll,phi_1), # this node phase acceptance [rad]
                omgl_0          = omgl_0,       # this node synchrotron oscillation [Hz]
                wmax            = wmax,         # this node w max on separatrix [1] (large amp. oscillations)
                zmax            = conv.DphiToz(-phisoll) # this node z max on separatrix [m] (large amp. oscillations -- Wrangler's approximation (pp.178) is good up to -58deg)
                )
        return res
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()
        fifo_z = kwargs['fifo_z']
        sfifo_z= kwargs['sfifo_z']
        fifo_xy= kwargs['fifo_xy']
        sfifo_xy=kwargs['sfifo_xy']

        rf_gap = self.gap_object      # this RF gap to use: can be self or others like OXAL_G or TTF_G

        lost=False
        tkin=rf_gap.particle.tkin

        res = self.waccept()
        dummy,phi_2,phisoll,phi_1 = res['phaseacc']   # dummy kept for old track_node
        conv = UTIL.WConverter(tkin,rf_gap.freq)
        Dphi=conv.zToDphi(new_point[Ktp.z])
        phi = phisoll+Dphi

        # longitudinal acceptance
        if not (phi_2 < phi and phi < phi_1):  # Wrangler's approximation (pp.178) is good up to -58deg
            fifo_z.append(f'loss (z) {new_point[Ktp.z]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        elif abs(new_point[Ktp.zp]) > res['Dp2pmax']:
            fifo_z.append(f'loss (zp) {new_point[Ktp.zp]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        # transverse apertures
        elif rf_gap.aperture != None and not (abs(new_point[Ktp.x]) < rf_gap.aperture or abs(new_point[Ktp.y]) < rf_gap.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost
class RFG(Node):   
    """  RF-gap of zero length for different kick gap-models """
    def __init__(self,label):
        super().__init__()
        self.aperture     = None
        self.cavlen       = None
        self.deltaW       = None
        self.dWf          = UTIL.FLAGS['dWf']
        self.EzPeak       = None
        self.freq         = None
        self.gap          = None
        self.label        = label
        self.length       = 0
        self.mapper       = None
        self.mapping      = UTIL.FLAGS.get('mapping')
        # self.matrix       = None
        self.particle     = UTIL.Proton(UTIL.PARAMS['injection_energy'])    
        self.particlef    = None
        self.position     = (0,0,0)
        self.phisoll      = None
        self.SFdata       = None
        self.sec          = None
        self.ttf          = None
        self.viseo        = 0.25

    def register(self,mapper):
        self.mapper  = mapper
        self.mapper.register(self)
        pass
    def configure(self,**kwargs): 
        if self.mapping in ['t3d','oxal','base','ttf']:
            self.aperture  = kwargs.get('aperture')
            self.cavlen    = kwargs.get('cavlen')          # [m] cavity length
            self.EzPeak    = kwargs.get('EzPeak')*self.dWf # [MV/m]
            self.freq      = kwargs.get('freq')            # [Hz]  RF frequenz
            self.gap       = kwargs.get('gap')             # [m] rf-gap
            self.phisoll   = kwargs.get('phisoll')         # [radians] soll phase
            self.sec       = kwargs.get('sec')
            self.SFdata    = kwargs.get('SFdata')

            self.mapper.configure(**kwargs)
            pass
        elif mapping in ['simple','dyn']:
            raise(UserWarning(wrapRED(f'mapping not ready {mapping}')))
            sys.exit()
        else:
            raise(UserWarning(wrapRED(f'missing implementation {mapping}')))
            sys.exit()
    @property
    def isAccelerating(self):
        return self.mapper.isAccelerating()
    @property
    def toString(self):
        return self.mapper.toString()
    def adjust_energy(self, tkin):
        self.mapper.adjust_energy(tkin)
        return self
    def map(self,i_track):
        return self.mapper.map(i_track)
    def make_slices(self, anz=1):
        return [self]
    def waccept(self):
        return self.mapper.waccept()
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()            # calling bunch.Tpoint() object

        fifo_z   = kwargs['fifo_z']   # collect losses in tracker.Fifo objects
        sfifo_z  = kwargs['sfifo_z']
        fifo_xy  = kwargs['fifo_xy']
        sfifo_xy = kwargs['sfifo_xy']

        lost   = False
        mapper = self.mapper
        tkin   = mapper.particle.tkin

        res  = self.waccept()    # test energy acceptance
        (dummy,phi_2,phisoll,phi_1) = res['phaseacc']   # dummy kept for old track_node
        conv = UTIL.WConverter(tkin, mapper.freq)
        Dphi = conv.zToDphi(new_point[Ktp.z])
        phi  = mapper.phisoll+Dphi

        # longitudinal acceptance
        if not (phi_2 < phi and phi < phi_1):  # Wrangler's approximation (pp.178) is good up to -58 [deg]
            fifo_z.append(f'loss (z) {new_point[Ktp.z]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        elif abs(new_point[Ktp.zp]) > res['Dp2pmax']:
            fifo_z.append(f'loss (zp) {new_point[Ktp.zp]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        # transverse apertures
        elif mapper.aperture != None and not (abs(new_point[Ktp.x]) < mapper.aperture or abs(new_point[Ktp.y]) < mapper.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost
class RFC(Node):   
    """ Rf cavity as product DKD*RFG*DKD 
    # def __init__(self, label, EzPeak, phisoll, gap, cavlen,freq, SFdata=0, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=UTIL.FLAGS['dWf'], mapping='t3d'):
    def __init__(self, label,mapping):
        super().__init__()

        self.viseo        = 0.25
        self.label        = label
        self.mapping      = mapping
        self.length       = 0.
        self.accelerating = True
        self.particle     = None
        self.mapper       = None
        self.ttf          = None
        self.deltaw       = None
        self.particlef    = None
        self.matrix       = None

        # self.dr = DKD(label="D-D",particle=particle,position=position,length=self.cavlen/2.,aperture=aperture)

    def register_mapping(self,mapper):
        self.mapper = mapper
        mapper.accept_register(self)

    def configure(self,**kwargs): 
        self.particle = kwargs['particle']  
        if self.mapping == "t3d":
            self.mapper.configure(**kwargs)

        elif self.mapping == 'simple':
            raise(UserWarning(wrapRED('mising implementation')))
            sys.exit()

        elif self.mapping == "oxal":
            raise(UserWarning(wrapRED('mising implementation')))
        
        else:
            raise(UserWarning(wrapRED('mising implementation')))

    def adjust_energy(self, tkin):
        self.mapper.adjust_energy(tkin)
        at_exit        = self.mapper.values_at_exit()
        self.ttf       = at_exit['ttf']
        self.deltaw    = at_exit['deltaw']
        self.particlef = at_exit['particlef']
        self.matrix    = at_exit['matrix']
        return self
    def map(self,i_track):
        return self.mapper.map(i_track)
    def make_slices(self, anz=1):
        return [self]
    def waccept(self):
        return self.mapper.waccept() """

def K(gradient, particle):
    """ quad strength K[1/m**2] for protons, gradient[T/m] """
    return 0.31952 * gradient/particle.gamma_beta

class TestElementMethods(unittest.TestCase):
    def Matrix(self,v):
        matrix = NP.eye(MDIM,MDIM)
        matrix[0,0] = v[0]
        matrix[0,1] = v[1]
        matrix[1,0] = v[2]
        matrix[1,1] = v[3]
        matrix[2,2] = v[4]
        matrix[2,3] = v[5]
        matrix[3,2] = v[6]
        matrix[3,3] = v[7]
        matrix[4,4] = v[8]
        matrix[4,5] = v[9]
        matrix[5,4] = v[10]
        matrix[5,5] = v[11]
        return matrix
    def test_Node(self):
        """ testing the * operator for NODE objects"""
        print("\b----------------------------------------test_Node")
        a = self.Matrix((1,2,3,4,1,2,4,5,1,2,5,6))
        b = self.Matrix((1,0,0,1,1,0,0,1,1,0,0,1))
        b = 2*b
        c = self.Matrix((1,1,0,1,1,1,0,1,1,1,0,1))
        A = Node()
        A.matrix = a
        A.label = 'Node-A'
        A.length = 1.
        # print(A.toString())

        B = Node()
        B.matrix = b
        B.label = 'Node-B'
        B.length = 1.
        # print(B.toString())
        AB = A*B
        # print(AB.toString())
        BA = B*A
        # print(BA.toString())

        self.assertTrue(NP.array_equal(AB.matrix,a.dot(b)),msg='AB * ab')
        self.assertTrue(NP.array_equal(BA.matrix,b.dot(a)),msg='BA * ba')
        self.assertTrue(NP.array_equal(AB.matrix,BA.matrix),msg='A * B')

        C = Node()
        C.matrix = c
        C.label = 'Node-C'
        C.length = 1.
        # print(C.toString())
        AC = A*C
        # print(AC.toString())
        CA = C*A
        # print(CA.toString())

        self.assertFalse(NP.array_equal(AC.matrix,CA.matrix))
    def test_QF_Node_make_slices(self):
        print("\b----------------------------------------test_QF_Node_make_slices")
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QFnode = QF("QFoc", gradient, particle=p, length=length)
        # slices = QFnode.make_slices(anz=anz)
        anz=5
        self.assertEqual(len(QFnode.make_slices(anz=anz)),anz)
        for i in range(anz):
            self.assertEqual(QFnode.make_slices(anz=anz)[i].length,length/anz)
        anz=1
        self.assertEqual(len(QFnode.make_slices(anz=anz)),anz)
        for i in range(anz):
            self.assertEqual(QFnode.make_slices(anz=anz)[i].length,length/anz)
        anz=0
        self.assertEqual(len(QFnode.make_slices(anz=anz)),1)
        self.assertEqual(QFnode.make_slices(anz=anz)[0].length,length)
    def test_I_Node(self):
        print("\b----------------------------------------test_I_Node")
        i_matrix = NP.eye(MDIM,MDIM)
        Inode    = I("IN1")
        II       = Inode*Inode
        self.assertEqual(Inode.type,"I",'type')
        self.assertEqual(Inode.label,"IN1","label")
        self.assertEqual(Inode.length,0.,"length")
        self.assertTrue(NP.array_equal(Inode.matrix,i_matrix),"matrix")
        self.assertTrue(NP.array_equal(II.matrix,i_matrix.dot(i_matrix)),"I*I")
        self.assertEqual(II.label,"IN1*IN1")
        # NOTE: type( Child * Child ) = 'Node'
        self.assertTrue(II.type == "Node")
    def test_D_Node(self):
        print("\b----------------------------------------test_D_Node")
        l = 10.
        p = UTIL.Proton(50.)
        Dnode = D("Drift",length=l,aperture=5.)
        self.assertEqual(Dnode.length,10.)
        g = p.gamma
        d_matrix = NP.eye(MDIM,MDIM)
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,DSKOO] = l
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"matrix")
        tkin = 200.
        Dnode = Dnode.adjust_energy(tkin)
        p = p(tkin)
        g = p.gamma
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"adjust_energy")
        l = l/2.
        Dnode = Dnode.shorten(l)
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,DSKOO] = l
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"shorten")
        self.assertTrue(Dnode.label == "Drift")
        self.assertTrue(Dnode.type == "D")
        self.assertEqual(Dnode.particle.tkin,p.tkin,"tkin")
    def test_Particle_and_K(self):
        print("\b----------------------------------------test_Particle_and_K")
        p = UTIL.Proton(50.)
        gradient =3.
        K0 = K(gradient,p)
        self.assertAlmostEqual(p.tkin, 50., delta=1e-3)
        self.assertAlmostEqual(p.gamma, 1.0533, delta=1e-4)
        self.assertAlmostEqual(p.beta, 3.1405e-1, delta=1e-4)
        self.assertAlmostEqual(p.gamma_beta, 3.3079e-1, delta=1e-4)
        self.assertAlmostEqual(p.brho, 1.0353, delta=1e-4)
        self.assertAlmostEqual(K0, 2.8978, delta=1e-4)
    def test_QF_Node(self): 
        print("\b----------------------------------------test_QF_Node")
        def matrix(gradient,length,particle,thin):
            x=0; xp=1; y=2; yp=3; z=4; zp=5; T=6; dT=7; S=8; dS=9
            k02  = K(gradient,particle)
            k0w  = M.sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[x,x]  = M.cos(dphi);        mx[x,xp]  = M.sin(dphi)/k0w       
                mx[xp,x] = -k0w*M.sin(dphi);   mx[xp,xp] = mx[x,x]
                mx[y,y]  = M.cosh(dphi);       mx[y,yp]  = M.sinh(dphi)/k0w        
                mx[yp,y] = +k0w*M.sinh(dphi);  mx[yp,yp] = mx[y,y]
            else:
                l = length
                f = 1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QFnode = QF("QFoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
        self.assertEqual(QFnode.viseo,+0.5,"viseo")
        self.assertTrue(QFnode.matrix[XPKOO,XKOO] < 0.)
        self.assertTrue(QFnode.matrix[YPKOO,YKOO] > 0.)

        p = UTIL.Proton(100.)
        QFnode = QFnode.adjust_energy(100.)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.adjust_energy")

        length = 0.5
        QFnode = QFnode.shorten(length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.shorten")

        gradient = 1.0; length = 0.15; p = UTIL.Proton(100.)
        QFnode = QF("QFfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
    def test_QD_Node(self): 
        print("\b----------------------------------------test_QD_Node")
        def matrix(gradient,length,particle,thin):
            x=0; xp=1; y=2; yp=3; z=4; zp=5; T=6; dT=7; S=8; dS=9
            k02  = K(gradient,particle)
            k0w  = M.sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[y,y]  = M.cos(dphi);        mx[y,yp]  = M.sin(dphi)/k0w        
                mx[yp,y] = -k0w*M.sin(dphi);   mx[yp,yp] = mx[y,y]
                mx[x,x]  = M.cosh(dphi);       mx[x,xp]  = M.sinh(dphi)/k0w        
                mx[xp,x] = +k0w*M.sinh(dphi);  mx[xp,xp] = mx[x,x]
            else:
                l = length
                f = -1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QDnode = QD("QDfoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
        self.assertTrue(QDnode.matrix[XPKOO,XKOO] > 0.)
        self.assertTrue(QDnode.matrix[YPKOO,YKOO] < 0.)

        length = 0.5
        QDnode = QDnode.shorten(length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"QD.shorten")

        self.assertEqual(QDnode.label,'QDfoc',"__getitem__")
        self.assertEqual(QDnode.type,'QD',"__getitem__")
        self.assertEqual(QDnode.viseo,-0.5,"viseo")
        QDnode.viseo = 0.75
        self.assertEqual(QDnode.viseo,0.75,"__setitem__")

        gradient = 1.0; length = 0.15; p = UTIL.Proton(100.)
        QDnode = QD("QDfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
    def test_Wille_Nodes(self):
        print("\b----------------------------------------test_Wille_Nodes")
        kqf  = -1.20    # [1/m**2]
        kqd  = -kqf
        lqf  = 0.20       # [m]
        lqd  = 0.40       # [m]
        rhob = 3.8197     # [m]
        lb   = 1.50       #[m]
        phib = 11.25      # [deg]
        ld   = 0.55       # [m]
        p    = UTIL.Proton(100.)
        gradf = kqf*p.brho
        gradd = kqd*p.brho

        """ Wille's Fokussierende Quadrupole = MQF """
        MQF = QF("MQF",gradf,particle=p,length=lqf)
        mqf = NP.eye(MDIM,MDIM)
        mqf[XKOO,XKOO]  = 0.9761;     mqf[XKOO,XPKOO]  = 0.1984
        mqf[XPKOO,XKOO] = -0.2381;    mqf[XPKOO,XPKOO] = mqf[XKOO,XKOO]
        mqf[YKOO,YKOO]  = 1.0241;     mqf[YKOO,YPKOO]  = 0.2016
        mqf[YPKOO,YKOO] = 0.2419;     mqf[YPKOO,YPKOO] = mqf[YKOO,YKOO]
        mqf[ZKOO,ZPKOO] = lqf/p.gamma**2
        mqf[EKOO,EKOO]  = 1.;         mqf[EKOO,DEKOO]  = 0.     # tkin
        mqf[SKOO,SKOO]  = 1.;         mqf[SKOO,DSKOO]   = lqf     # s
        # print("\nmqf");print(mqf);print("QF.matrix");print(MQF.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mqf[i,j],MQF.matrix[i,j],msg="MQF",delta=1e-4)

        """ Wille's Defokussierende Quadrupole = MQD """
        MQD = QD("MQD",gradd,particle=p,length=lqd)
        mqd = NP.eye(MDIM,MDIM)
        mqd[XKOO,XKOO]  = 1.0975;     mqd[XKOO,XPKOO]  = 0.4129
        mqd[XPKOO,XKOO] = 0.4955;     mqd[XPKOO,XPKOO] = mqd[XKOO,XKOO]
        mqd[YKOO,YKOO]  = 0.9055;     mqd[YKOO,YPKOO]  = 0.3873
        mqd[YPKOO,YKOO] = -0.4648;    mqd[YPKOO,YPKOO] = mqd[YKOO,YKOO]
        mqd[ZKOO,ZPKOO] = lqd/p.gamma**2
        mqd[EKOO,EKOO]  = 1.;         mqd[EKOO,DEKOO]  = 0.     # tkin
        mqd[SKOO,SKOO]  = 1.;         mqd[SKOO,DSKOO]   = lqd     # s
        # print("\nmqd");print(mqd);print("MQD.matrix");print(MQD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mqd[i,j],MQD.matrix[i,j],msg="MQD",delta=1e-4)

        """ Wille's Dipole = MB """
        MB = SD("MB",2*phib,rhob,particle=p)
        mb = NP.eye(MDIM,MDIM)
        mb[XKOO,XKOO]  = 0.9239;      mb[XKOO,XPKOO]   = 1.4617;        mb[XKOO,ZPKOO] = 0.2908
        mb[XPKOO,XKOO] = -0.1002;     mb[XPKOO,XPKOO]  = mb[XKOO,XKOO]; mb[XPKOO,ZPKOO] = 0.3827
        mb[YKOO,YKOO]  = 1.0;         mb[YKOO,YPKOO]   = 1.5
        mb[YPKOO,YKOO] = -0.0;        mb[YPKOO,YPKOO]  = mb[YKOO,YKOO]
        mb[ZKOO,XKOO]  = -0.3827;     mb[ZKOO,XPKOO]   = -0.2908;       mb[ZKOO,ZPKOO] = 1.1867
        mb[EKOO,EKOO]  = 1.;          mb[EKOO,DEKOO]   = 0.     # tkin
        mb[SKOO,SKOO]  = 1.;          mb[SKOO,DSKOO]    = lb     # s
        # print("\nmb");print(mb);print("MB.matrix");print(MB.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mb[i,j],MB.matrix[i,j],msg="MB",delta=1e-4)

        """ Wille's Kantenfokusierung = MEB """
        MEB = Wedge(phib,rhob, t3d_wedge=False)
        meb = NP.eye(MDIM,MDIM)
        meb[XPKOO,XKOO] = 0.0521; meb[YPKOO,YKOO] = -meb[XPKOO,XKOO]
        # print("\nmeb");print(meb);print("MEB.matrix");print(MEB.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(meb[i,j],MEB.matrix[i,j],msg="MEB",delta=1e-2)

        """ Wille's Driftstrecken = MD """
        MD = D("MD",particle=p,length=ld)
        md = NP.eye(MDIM,MDIM)
        md[XKOO,XPKOO] = ld; md[YKOO,YPKOO] = md[XKOO,XPKOO]
        md[ZKOO,ZPKOO] = ld/p.gamma**2
        md[SKOO,DSKOO] = ld
        # print("\nmd");print(md);print("MD.matrix");print(MD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(md[i,j],MD.matrix[i,j],msg="MD",delta=1e-4)

        """ Wille's gesamte Zelle mit Sektordipol und Kantenfokussierung = MZ """
        MZ = MQF*MD*MEB*MB*MEB*MD*MQD*MD*MEB*MB*MEB*MD*MQF
        # MZ = MQF*MD*MEB*MB*MEB*MD*MQD   #*MD*MEB*MB*MEB*MD*MQF
        mz = NP.eye(MDIM,MDIM)
        mz[XKOO,XKOO]  = 0.0808;   mz[XKOO,XPKOO]  = 9.7855;        mz[XKOO,ZPKOO]  = 3.1424
        mz[XPKOO,XKOO] = -0.1015;  mz[XPKOO,XPKOO] = mz[XKOO,XKOO]; mz[XPKOO,ZPKOO] = 0.3471
        mz[YKOO,YKOO]  = -0.4114;  mz[YKOO,YPKOO]  = 1.1280
        mz[YPKOO,YKOO] = -0.7365;  mz[YPKOO,YPKOO] = mz[YKOO,YKOO]
        mz[ZKOO,XKOO]  = -0.34707; mz[ZKOO,XPKOO]  = -3.1424;       mz[ZKOO,ZPKOO]  = 4.1844
        mz[SKOO,DSKOO]  = 6.0
        # print("\nmz");print(mz);print("MZ.matrix");print(MZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MZ.matrix[i,j],msg="MZ",delta=1e-4)
        """ Wille's gesamte Zelle mit RD-class """
        WEDGE = MEB
        MRD   = RD("RD",2*phib,rhob,WEDGE,particle=p)
        MRDZ  = MQF*MD*MRD*MD*MQD*MD*MRD*MD*MQF
        # print("\nmz");print(mz);print("MRDZ.matrix");print(MRDZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MRDZ.matrix[i,j],msg="MRDZ == MZ?",delta=1e-4)
    def test_SD_and_RD_Node_make_slices(self):
        print("\b----------------------------------------test_SD_and_RD_Node_make_slices")
        rhob   = 3.8197   # [m]
        phib   = 11.25    # [deg]
        p      = UTIL.Proton(100.)

        """ slice,recombine and adjust SD """
        sd     = SD("SD",2*phib,rhob,particle=p)
        slices = sd.make_slices(anz=3)
        # for item in slices: print(item.matrix); print()
        sdx = slices[0]         # combine the slices again
        for i in range(1,len(slices)):
            sdx = sdx * slices[i]
        # print("\nsd"); print(sd.matrix); print("\nsdx"); print(sdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(sd.matrix[i,j],sdx.matrix[i,j],msg='sd == sdx?',delta=1.e-4)

        sd_adjusted = sd.adjust_energy(tkin=p.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(sd.matrix[i,j],sd_adjusted.matrix[i,j],'sd == sd_adjusted?')

        """ slice,recombine and agjust RD """
        wedge =  Wedge(phib,rhob,t3d_wedge=True)
        rd = RD("RD",2*phib,rhob,wedge,particle=p)
        slices = rd.make_slices(anz=5)
        # for item in slices: print(item.matrix); print()
        rdx = slices[0]
        for i in range(1,len(slices)):
            rdx = rdx * slices[i]
        # print("\nrd"); print(rd.matrix); print("\nrdx"); print(rdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(rd.matrix[i,j],rdx.matrix[i,j],msg='rd == rdx?',delta=1.e-4)

        rd_adjusted = rd.adjust_energy(tkin=p.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(rd.matrix[i,j],rd_adjusted.matrix[i,j],'rd == rd_adjusted?')
    def test_SD_and_RD_Node_adjust_energy(self):
        print("\b----------------------------------------test_SD_and_RD_Node_adjust_energy")
        rhob   = 3.8197   # [m]
        phib   = 11.25    # [deg]
        p100      = UTIL.Proton(100.)
        p50       = UTIL.Proton(50.)

        """ adjust SD  100->50->100 """
        sd = SD("SD",2*phib,rhob,particle=p100)
        sd_adjusted = sd.adjust_energy(tkin=p50.tkin)
        sd_100 = sd_adjusted.adjust_energy(tkin=p100.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(sd.matrix[i,j],sd_100.matrix[i,j],'sd == sd_100?')

        """  adjust RD 100->50->100 """
        wedge =  Wedge(phib,rhob,t3d_wedge=True)
        rd = RD("RD",2*phib,rhob,wedge,particle=p100)
        rd_adjusted = rd.adjust_energy(tkin=p50.tkin)
        rd_100 = rd_adjusted.adjust_energy(tkin=p100.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(rd.matrix[i,j],rd_100.matrix[i,j],'rd == rd_100?')
    def test_GAP_Node(self):
        print("\b----------------------------------------test_GAP_Node")
        EzAvg   = 2.1             #[MV/m]
        phisoll = M.radians(-30.)
        gap     = 0.022           #[m]
        freq    = 816.e6          #[Hz]
        gap = GAP("GAP",EzAvg,phisoll,gap,freq,dWf=1) # tkin=500 default
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.03766,delta=1.e-4)
        gap = gap.adjust_energy(6.) # tkin=6
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.023817,delta=1.e-4)
    def test_RFG_Node_with_T3D_G_map(self):
        print("\b----------------------------------------test_RFG_Node_with_T3D_G_map")
        ID        = 'RFG mit T3D_G'
        EzPeak    = 2.1
        phisoll   = M.radians(-30.)
        gap_parameters = dict(
            EzPeak    = EzPeak,
            phisoll   = phisoll,         # [radians] requested soll phase
            gap       = 0.044,
            cavlen    = 0.064,
            freq      = 816.e6,            # [Hz]  requested RF frequenz
            particle  = Proton(50.),
            position  = (0,0,0),
            aperture  = 0.011,
            sec       = 'xx'
        )
        instance = RFG(ID)
        instance.register_mapping(T3D_G())
        instance.configure(**gap_parameters)
        
        self.assertEqual(instance.mapping,'t3d')
        self.assertEqual(instance.label,'RFG mit T3D_G')
        self.assertEqual(instance.mapper.label,'T3D_G')
        self.assertEqual(instance.mapper.EzPeak,EzPeak)
        self.assertEqual(instance.mapper.phisoll,M.radians(-30.))
        self.assertAlmostEqual(instance.mapper.deltaW,0.062206,delta=1.e-4)
        self.assertAlmostEqual(instance.mapper.matrix[EKOO,DEKOO],0.062206,delta=1.e-4)
        self.assertEqual(instance.length,0.)
        self.assertEqual(instance.mapper.matrix[SKOO,DSKOO],0.)
        instance.adjust_energy(250.)
        self.assertAlmostEqual(instance.mapper.deltaW,0.0751,delta=1.e-4)
        self.assertAlmostEqual(instance.matrix[EKOO,DEKOO],0.0751,delta=1.e-4)

if __name__ == '__main__':
    UTIL.FLAGS['verbose'] = 3
    unittest.main()
