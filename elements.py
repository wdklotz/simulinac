#!/Users/klotz/anaconda3/bin/python3.6
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
from math import sqrt, sinh, cosh, sin, cos, tan, modf, pi, radians, ceil
from copy import copy, deepcopy
import numpy as NP
import pprint, inspect
import unittest

from setutil import PARAMS, FLAGS, Particle
from setutil import WConverter, dictprnt, objprnt, Proton, Electron
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO, MDIM
from setutil import dBdxprot, scalek0prot, k0prot, I0, I1, arrprnt, Ktp
from Ez0     import SFdata
from TTFG    import _TTF_G
from DynacG  import _DYN_G
from OXAL    import _OXAL

def PRINT_PRETTY(obj=None):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    if obj != None: pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj=None):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

twopi = 2.*pi     # used about everywhere
# numpy pretty printing
NP.set_printoptions(linewidth = 132, formatter = {'float': '{:>8.5g}'.format})

class OutOfRadialBoundEx(Exception):
    def __init__(self,max_r,ID=''):
        self.message = "OutOfRadialBoundEx: in '{}' out of {} [cm] max radial excursion.".format(ID,max_r*100.)
#------- The mother of all lattice elements (a.k.a. matrices)
#------- The mother of all lattice elements (a.k.a. matrices)
#------- The mother of all lattice elements (a.k.a. matrices)
class Node(object):
    """ Base class for transfer matrices (linear map)
        ii)  is a dictionary (DictObject base class)
        ii)  each instance holds its copy of the refrence particle (self.particle)
    """
    slice_min = 0.001                   # default - minimal slice length
    def __init__(self):
        self.type      = self.__class__.__name__      # self's node type
        self.particle  = None      # !!!IMPORTANT!!! local copy of the particle object
        self.matrix    = None      # MDIMxMDIM zero matrix used here
        self.position  = None      # [entrance, middle, exit]
        self.length    = None      # default - thin
        self.label     = None      # default - unlabeled
        self.aperture  = None      # default - infinite aperture
        self.next      = None      # right link
        self.prev      = None      # left link
        self.viseo     = None      # default - invisible
    def __getitem__(self,k):
        return self.__dict__[k]
    def __setitem__(self,k,v):
        self.__dict__[k] = v
    def toString(self):
        ret = repr(self)
        for k,v in self.__dict__.items():
            ret+= '\n{}:{}'.format(k,v)
        return ret
    def adjust_energy(self, tkin):
        """ nothing to adjust """
        return self
    @property
    def twiss(self):
        return self['twiss']
    @property
    # def particlef(self):
    #      # !!!IMPORTANT!!! return a copy with updated energy
    #     return copy(self.particle)(self.particle.tkin + self.deltaW)
    def __call__(self, n = MDIM, m = MDIM):
        # return upper left n, m submatrix
        return self.matrix[:n, :m]
    def __mul__(self, other):
        """ define the (*) operator for Node objects """
        res = Node()
        if (self.label == ''):
            res.label = other.label
        else:
            res.label = self.label+'*'+other.label
        res.length = self.length + other.length
        """ matrix product """
        res.matrix = NP.dot(self.matrix, other.matrix)
        return res
    def prmatrix(self):
        n  = 1000
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
    def reverse(self):
        raise RuntimeError('Node:reverse not implemented!')
        sys.exit()
    def inverse(self):
        raise RuntimeError('Node:inverse not implemented!')
        sys.exit(1)
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
    def make_slices(self, anz = PARAMS['nbslices']):
        # ignore the very small rest
        mr = None
        slices = []

        if self.length == 0.:
            # nothing todo for zero length element
            slices.append(self)

        else:
            step = self.length/anz      # calc step size
            if step < self['slice_min']:
                step  = self['slice_min']

            (step_fraction_part, step_int_part) = modf(self.length/step)

            rest = step * step_fraction_part
            # shorten element to step length
            mx   = self.shorten(step)
            if rest > 1.e-3:
                mr = self.shorten(rest)
            elif rest < 0.:
                raise RuntimeError('FATAL: negative resting step size when stepping through - STOP')
                sys.exit(1)

            for i in range(int(step_int_part)):
                slices.append(mx)
        if mr != None: slices += [mr]
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
    def sigma_beam(self, steps=1, sg=None):
        """ 
        Track the Sigma object through a node
            twiss vector: twv  = NP.array([betax,alphax,gammax,b..y,a..y,g..y,b..z,a..z,g..z])
            *) input: sg = SIGMA(twv,epsx,epsy,epsz) Sigma object
        """
        si,sm,sf = self.position # entrance
        sigmas   = []
        s        = si
        slices = self.make_slices(anz = steps)
        for slice in slices:
            # next_SIGMA = R * SIGMA * transpose(R)
            sgf = sg.RSRT(slice)
            # emmitance grow ?
            if isinstance(slice,RFG) and FLAGS['egf']: # loop slices
                sgf = sgf.apply_eg_corr(rf_gap=slice, sigma_i=sg, delta_phi=PARAMS['Dphi0'])
            s += slice.length
            sigmas.append((sgf,s))
            sg = sgf

        # averages
        av = []
        for sig,s in sigmas:
            v = sig.sigv().tolist()
            av.append(v)
        avarr = NP.array(av)
        avm = NP.mean(avarr,axis=0)
        # sgxm, sgym are mean values of sigmas over slices
        sgxm  = avm[0]; sgxpm = avm[1]
        sgym  = avm[2]; sgypm = avm[3]
        sigxy = (sgxm,sgxpm,sgym,sgypm)
        # each node has its tuple of average sigmas     
        self['sigxy'] = sigxy
        return sigmas
    def map(self, i_track):
        """ Linear mapping of trajectory from (i) to (f) """
        f_track = copy(i_track)
        f_track = NP.dot(self.matrix,f_track)

        # for DEBUGGING
        if 0:
            f = f_track.copy() # !!!IMPORTANT!!!
            for i in range(len(f_track)-4):
                f[i]  = f[i]*1.e3
            arrprnt(f, fmt = '{:6.3g},', txt = 'matrix_map: ')
        return f_track
    def soll_map(self, i_track):
        f_track = copy(i_track)
        f_track[EKOO] += i_track[DEKOO]*self.matrix[EKOO,DEKOO]
        return f_track
    def shorten(self, length):
        """ nothing to shorten """
        return self
class I(Node):
    """  Unity matrix: the unity Node """
    def __init__(self, label):
        super().__init__()
        self.label    = label
        self.matrix   = NP.eye(MDIM,MDIM) # set the NODE's member variable
        self.length   = 0.
class MRK(I):
    """ Marker node (a.k.a element): owns a list of agents that do the actions """
    def __init__(self, label, agents=[], particle=PARAMS['sollteilchen'], position=[0, 0, 0]):
        super().__init__(label)
        self.agents = agents   # the agent-list
        self.particle = copy(particle)
        self.position = position
        self['viseo'] = 4.
    def add(self,agent):
        self.agents.append(agent)
    def do_actions(self):
        """ invoke all actions bound to this marker """
        for agent in self.agents:
            agent.do_action()
    def adjust_energy(self, tkin):
        adjusted = MRK(self.label, agents=self.agents, particle=self.particle(tkin), position=self.position)
        return adjusted
class D(Node):
    """ 
    Trace3D drift space 
    """
    def __init__(self, label, particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0.,aperture=PARAMS['aperture']):
        super().__init__()
        self.label    = label
        self.particle = copy(particle)
        self.position = position
        self.length   = length
        self.aperture = aperture
        self.matrix   = NP.eye(MDIM,MDIM)
        m = self.matrix 
        g = self.particle.gamma
        m[XKOO, XPKOO] = m[YKOO, YPKOO] = self.length
        m[ZKOO, ZPKOO] = self.length/(g*g)
        m[SKOO, LKOO]  = self.length # Ds the longitudinal length increase

    def adjust_energy(self, tkin):
        adjusted = D(self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortend =  D(self.label, particle=self.particle, position=[0,0,0], length=length, aperture=self.aperture)
        return shortend
class QF(Node):
    """ 
    Trace3D focussing quad  !!! NEUES API !!!!
    """
    def __init__(self, label, grad, particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture']):
        super().__init__()
        self.label    = label
        self.grad     = abs(grad)                        # [T/m]
        self.particle = copy(particle)
        self.k02      = K(self.grad,self.particle)  # [1/m**2]
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
        k0w   = sqrt(self.k02)
        l = self.length
        k02 = self.k02
        phi = l*k0w
        f = 1./(k02*l)
        if f/l > 50.: self.thin = True
        DEBUG_OFF(F"{self.type}: k02={k02:.2f} \tf={f:.2f} >> l={l:.2f}? {self.thin}")
        if self.thin != True:
            """ thick quad """
            # focusing
            cf   = cos(phi)
            sf   = sin(phi)/k0w
            cfp  = -k0w*sin(phi)
            sfp  = cf 
            # defocusing
            cd  =  cosh(phi)
            sd = sinh(phi)/k0w
            cdp = k0w*sinh(phi)
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
        m[SKOO, LKOO]  = self.length # length increase
        return m
    def isThin(self):
        return self.thin
    def adjust_energy(self, tkin):
        adjusted = QF(self.label, self.grad, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortened = QF(self.label, self.grad, particle=self.particle, position=self.position, length=length, aperture=self.aperture)
        return shortened
class QD(QF):
    """ 
    Trace3D defocussing quad  !!! NEUES API !!!! 
    """
    def __init__(self, label, grad, particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture']):
        super().__init__(label, grad, particle=particle, position=position, length=length, aperture=aperture)
        self['viseo'] = -0.5

    def shorten(self, length):
        shortened = QD(self.label, self.grad, particle=self.particle, position=self.position, length=length,  aperture=self.aperture)
        return shortened
class SD(Node):
    """ Trace3d horizontal sector magnet. n=0 pure dipole. """
    def __init__(self, label, alpha, rho, n=0, particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=PARAMS['aperture']):
        super().__init__()
        self.label    = label
        self.alpha    = alpha   # [deg]
        self.rho      = rho     # [m]
        self.n        = n
        self.particle = copy(particle)
        self.position = position
        self.length   = None
        self.aperture = aperture
        self.viseo    = 0.25
        self.matrix   = self._mx()
    def _mx(self):
        m = NP.eye(MDIM,MDIM)
        beta  = self.particle.beta
        gamma = self.particle.gamma
        alpha = radians(self.alpha)      # [rad]
        rho   = self.rho                 # [m]
        h = alpha/(abs(rho*alpha))       # [1/m]
        kx = sqrt((1-self.n)*h**2)       # [1/m]
        ky = sqrt(self.n*h**2)
        Ds = abs(self.rho)*self.alpha         # [m]
        self.length = Ds
        cx = cos(kx*Ds)
        sx = sin(kx*Ds)
        cy = cos(ky*Ds)
        sy = sin(ky*Ds)

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

        m[SKOO, 
        LKOO]  = self.length  # length increase
        return m
    def adjust_energy(self, tkin):
        adjusted = SD(self.label,self.alpha,self.rho,self.n,particle=self.particle(tkin),position=self.position,aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        alpha = length/self.rho
        shortSD = SD(self.label,alpha,self.rho,self.n, particle=self.particle,position=self.position, aperture=self.aperture)
        return shortSD
    def make_slices(self, anz=PARAMS['nbslices']):
        shortSD = self.shorten(self.length/anz)
        slices = []          # wedge @ entrance
        for i in range(anz):
            slices.append(shortSD)
        DEBUG_OFF('slices {}'.format(slices))
        return slices
class RD(SD):
    """ Trace3D rectangular dipole x-plane """
    def __init__(self, label, alpha, rho, wedge, particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=PARAMS['aperture']):
        super().__init__(label, alpha, rho, particle=particle, position=position, aperture=aperture)
        self.wedge = wedge
        self.matrix = NP.dot(self.wedge.matrix,NP.dot(self.matrix,self.wedge.matrix))
    def make_slices(self, anz=PARAMS['nbslices']):
        slices = [self.wedge]          # wedge @ entrance
        shortSD = super().shorten(self.length/anz)
        for i in range(anz):
            slices.append(shortSD)
        slices.append(self.wedge)
        DEBUG_OFF('slices {}'.format(slices))
        return slices
class _wedge(Node):
    """ 
    Trace3d wedge Rxx and Ryy simplified
    """
    def __init__(self, beta, rho, t3d_wedge=True):
        super().__init__()
        self.label = "W"
        self.length = 0.
        g = 0.050    # fixed gap of 50 mm assumed
        K1=0.45
        K2=2.8
        psi  = K1*g/rho*((1+sin(beta)**2)/cos(beta))*(1-K1*K2*(g/rho)*tan(beta)) if t3d_wedge else 0.
        mxpx = tan(beta)/rho
        mypy = tan(beta-psi)/rho
        mx   = NP.eye(MDIM,MDIM)
        mx[XPKOO, XKOO] = mxpx
        mx[YPKOO, YKOO] = -mypy
        self.matrix = mx
class GAP(Node):
    """ Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    ... nicht sehr nuetzlich: produziert keine long. Dynamik wie Trace3D RFG!  """
    def __init__(self, label, EzAvg, phisoll, gap, freq, particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=None, dWf=FLAGS['dWf']):
        """ EzAvg [MV/m], phisoll [rad], gap [m], freq [Hz] """
        super().__init__()
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
        lamb   = PARAMS['clight']/self.freq            # [m] wellenlaenge
        beta   = self.particle.beta                    # beta Einstein
        ttf    = self.ttf(beta,lamb,self.gap)          # time-transition factor
        bg     = self.particle.gamma_beta
        m0c2   = self.particle.e0
        deltaW = E0L*ttf*cos(self.phisoll)                   # delta-W T.Wrangler pp.221
        m      = NP.eye(MDIM,MDIM)
        cyp = cxp = -pi*E0L*ttf*sin(self.phisoll)/(m0c2*lamb*bg**3)
        m[XPKOO, XKOO] = cxp
        m[YPKOO, YKOO] = cyp
        m[EKOO, DEKOO] = deltaW    # energy increase
        return deltaW,m
    def ttf(self, beta, lamb, gap):
        """ ttf-factor nach Panofsky (Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
        x = pi*gap/(beta*lamb)
        ttf = NP.sinc(x/pi)
        return ttf
    def adjust_energy(self, tkin):
        adjusted = GAP(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,particle=self.particle(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted
class RFG(Node):
    """  Wrapper to zero length RF kick gap-models """
    def __init__(self, label, EzAvg, phisoll, gap, freq, mapping='t3d', SFdata = None, particle=PARAMS['sollteilchen'], position=[0,0,0], aperture=None, dWf=FLAGS['dWf']):
        super().__init__()
        self.label     = label
        self.EzAvg     = EzAvg*dWf          # [MV/m] average gap field
        self.phisoll   = phisoll            # [radians] soll phase
        self.freq      = freq               # [Hz]  RF frequenz
        self.omega     = twopi*self.freq
        self.lamb      = PARAMS['clight']/self.freq
        self.gap       = gap                # [m] rf-gap
        self.mapping   = mapping            # map model
        self.SFdata    = SFdata             # SuperFish data
        self.particle  = copy(particle)
        self.position  = position
        self.aperture  = aperture
        self.length    = 0.
        self.dWf       = dWf                 # dWf=1 wirh acceleration else 0
        self.viseo     = 0.25
        self.matrix    = None
        self.ttf       = None
        """ RFG 'hosts' gap_model - which is a ref. to the 'client' class, i.e. _T3D_G as default """
        self.gap_model = _T3D_G(host=self) 
        """ set dispatching to gap models """
        if self.mapping == 't3d': pass
        elif self.mapping == 'simple' or self.mapping == 'base':
            self.gap_model = _PYO_G(host=self) # PyOrbit gap-models w/o SF-data
        elif self.mapping == 'ttf':
            self.gap_model = _TTF_G(host=self) # 3 point TTF-RF gap-model with SF-data  (A.Shishlo/J.Holmes)
        elif self.mapping == 'oxal':
            self.gap_model = _OXAL(host=self) # openXAL gap-model with SF-data  (A.Shishlo/J.Holmes)
        else:
            print(F"INFO: RFG is a kick-model and does not work with {self.mapping} mapping! Use one of [t3d,simple,base,ttf,oxal].")
    def adjust_energy(self, tkin):
        adjusted = RFG(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,mapping=self.mapping,SFdata=self.SFdata,particle=self.particle(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
        return adjusted

    """ delegate mapping to gap-model """
    def map(self, i_track):
        return self.gap_model.map(i_track)
    def soll_map(self, i_track):
        return self.gap_model.soll_map(i_track)
class _T3D_G(Node):   
    """ Mapping (i) to (f) in Trace3D zero length RF-Gap Model """
    def __init__(self, host=None):
        super().__init__()
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65) """
            x = pi*gap/(beta*lamb)
            return NP.sinc(x/pi)   # sinc(x) = sin(pi*x)/(pi*x)
        def mx(ttf, particlei, particlef, E0L, phisoll, lamb, deltaW,length):
            """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
            Wav     = particlei.tkin+deltaW/2.   # average tkin
            pav     = Proton(tkin=Wav)
            bav     = pav.beta
            gav     = pav.gamma
            m0c2    = pav.e0
            kz      = twopi*E0L*ttf*sin(phisoll)/(m0c2*bav*bav*lamb)
            ky      = kx = -0.5*kz/(gav*gav)
            bgi     = particlei.gamma_beta
            bgf     = particlef.gamma_beta
            bgi2bgf = bgi/bgf
            m       = NP.eye(MDIM,MDIM)
            m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
            m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
            m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
            # UPDATE NODE matrix with deltaW
            m[EKOO, DEKOO] = deltaW
            m[SKOO,LKOO]   = length
            return m
        # function body starts here -------------function body starts here -------------function body starts here -------------
        # function body starts here -------------function body starts here -------------function body starts here -------------
        # function body starts here -------------function body starts here -------------function body starts here -------------
        host.ttf       = ttf(host.lamb,host.gap,host.particle.beta)
        E0L            = host.EzAvg*host.gap
        host.deltaW    = E0L*host.ttf*cos(host.phisoll) # deltaW energy kick Trace3D
        tkin           = host.particle.tkin
        host.particlef = copy(host.particle)(tkin+host.deltaW) # !!!IMPORTANT!!! particle @ (f)
        host.matrix    = mx(host.ttf,host.particle,host.particlef,E0L,host.phisoll,host.lamb,host.deltaW,host.length)
    def map(self, i_track):
        """ Mapping from (i) to (f) with linear Trace3D matrix """
        f_track = copy(i_track)
        f_track = NP.dot(self.matrix,f_track)
        DEBUG_OFF('t3d-map {}'.format(f_track))
        return f_track
    def soll_map(self, i_track):
        si,sm,sf = self.position
        f_track = copy(i_track)
        f_track[EKOO] = self.particlef.tkin
        f_track[SKOO] = sf
        DEBUG_OFF('t3d-soll {}'.format(f_track))
        return f_track
class _PYO_G(object):
    """ 
    PyOrbit zero length RF gap-model (A.Shishlo,Jeff Holmes) 
    """
    def __init__(self, host=None):
        def ttf(lamb, gap, beta):
            """ Transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65) """
            x = gap/(beta*lamb)
            ttf = NP.sinc(x)   # sinc(x) = sin(pi*x)/(pi*x)
            return ttf
        # function body starts here -------------function body starts here -------------function body starts here -------------
        # function body starts here -------------function body starts here -------------function body starts here -------------
        # function body starts here -------------function body starts here -------------function body starts here -------------
        self.particle   = host.particle
        self.phis       = host.phis
        self.lamb       = host.lamb
        self.freq       = host.freq
        self.mapping    = host.mapping
        self.E0L        = host.EzAvg*host.gap
        self.ttf        = ttf(host.lamb, host.gap, host.particle.beta)
        self.qE0LT      = self.E0L*self.ttf
        # deltaW soll-energy increase Trace3D (same as Shishlo)
        self._deltaW    = self.qE0LT*cos(self.phis)
        self._particlef = copy(self.particle)(self.particle.tkin+self._deltaW)

        # UPDATE linear NODE matrix with deltaW
        host.matrix[EKOO, DEKOO] = self.deltaW

        if self.mapping == 'simple':
            self.which_map = self.simple_map
        elif self.mapping == 'base':
            self.which_map = self.base_map_1
            
    @property
    def deltaW(self):
        return self._deltaW
    @property
    def particlef(self):
        return self._particlef
    def map(self, i_track):
        return self.which_map(i_track)
    def soll_map(self, i_track):
        f_track       = copy(i_track)
        f_track[EKOO] = f_track[EKOO] + self.deltaW
        return f_track
    def simple_map(self, i_track):
        """ Mapping (i) to (f) in Simplified Matrix Model. (A.Shislo 4.1) """
        xi        = i_track[XKOO]       # [0]
        xpi       = i_track[XPKOO]      # [1]
        yi        = i_track[YKOO]       # [2]
        ypi       = i_track[YPKOO]      # [3]
        zi        = i_track[ZKOO]       # [4] z
        zpi       = i_track[ZPKOO]      # [5] Dp/p
        T         = i_track[EKOO]       # [6] kinetic energy SOLL
        S         = i_track[SKOO]       # [8] position SOLL

        particle  = self.particle
        m0c2      = particle.e0
        WIN       = particle.tkin
        betai     = particle.beta
        gammai    = particle.gamma
        gbi       = particle.gamma_beta
        deltaW    = self.deltaW
        lamb      = self.lamb
        qE0LT     = self.qE0LT
        phis      = self.phis
        

        DEBUG_OFF('simple_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,1.,phis))

        DW         = deltaW
        WOUT       = WIN + DW
        particlef  = copy(particle)(tkin = WOUT) # !!!IMPORTANT!!!
        betaf      = particlef.beta
        gammaf     = particlef.gamma
        gbf        = particlef.gamma_beta

        # the longitudinal 2x2 map (always linear!) A.Shishlo/J.Holmes (4.1.6-10)
        m11 = gbf/gbi
        m12 = 0.
        m21 = qE0LT*twopi/(lamb*betai)*sin(phis)
        m22 = 1.
 
        # Dp/p --> DW
        condPdT = m0c2*betai**2*gammai
        dwi     = condPdT*zpi  

        # long. (z,DW)
        zf  = m11*zi + m12*dwi
        dwf = m21*zi + m22*dwi

        # DW --> Dp/p
        condTdP = 1./(m0c2*betaf**2*gammaf)
        zfp     = dwf*condTdP

        # transverse (x',y')
        xpf  = gbi/gbf*xpi - xi * (pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf)) * sin(phis) # A.Shishlo/J.Holmes 4.1.11)
        ypf  = gbi/gbf*ypi - yi * (pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf)) * sin(phis)

        f_track = NP.array([xi, xpf, yi, ypf, zf, zfp, T, 1., S, 1.])

        # for DEBUGGING
        if 0:
            itr = i_track.copy()
            ftr = f_track.copy()
            for i in range(len(f_track)-4):
                itr[i]  = itr[i]*1.e3
                ftr[i]  = ftr[i]*1.e3
            arrprnt(itr, fmt = '{:6.3g},', txt = 'simple_map:i_track:')
            arrprnt(ftr, fmt = '{:6.3g},', txt = 'simple_map:f_track:')

        # the parent delegates reading these properties from here
        self._particlef = copy(particle)(particle.tkin + deltaW) # !!!IMPORTANT!!!
        return f_track
    def base_map_0(self, i_track):
        """alte map version bis 02.02.2022"""
        def DEBUG_TRACK(inout,track):
            print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        # function body ================= function body ================= function body ================= 
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
        frq      = self.freq
        lamb     = self.lamb
        phis     = self.phis
        qE0LT    = self.qE0LT
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)
        max_r  = 0.05              # max radial excursion [m]
        r      = sqrt(x**2+y**2)  # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(max_r,ID='_PYO_G:base_map')
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = I0(Kr)                               # bessel function I0
        i1     = I1(Kr)                               # bessel function I1
        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)
        # SOLL
        WIN       = tki                               # energy (i)
        DELTAW    = self.deltaW                       # energy kick
        WOUT      = WIN + DELTAW                      # energy (f) (4.1.6) A.Shishlo/J.Holmes
        # PARTICLE
        converter = WConverter(WIN,frq)
        # phin      = -z * twopi/(betai*lamb) + phis     # phase (i)  alte methode
        phin      = converter.zToDphi(z) + phis          # phase (i)
        deltaW    = qE0LT*i0*cos(phin)                   # energy kick
        # win     = (zp * (gammai+1.)/gammai +1.) * WIN  # energy (i) dp/p --> dT alte methode
        win       =  converter.Dp2pToW(zp) + WIN         # energy (i) dp/p --> dT
        wout      = win + deltaW                         # energy (f)   (4.2.3) A.Shishlo/J.Holmes
        dw        = wout - WOUT                          # d(deltaW)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particlef = copy(particle)(tkin = WOUT)       # !!!IMPORTANT!!! SOLL particle (f)
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
            xp  = gbi/gbf*xp - x/r*commonf*sin(phin)  # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbf*yp - y/r*commonf*sin(phin)
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x, xp, y, yp, z, zpf, T, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        if 0:
            arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
            arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        # the parent reads these attributes below
        self._particlef = particlef
        return f_track
    def base_map_1(self, i_track):
        """Neue map Version ab 03.02.2022 ist ein Remake um Korrecktheit der Rechnung zu testen."""
        def DEBUG_TRACK(inout,track):
            print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        # function body ================= function body ================= function body ================= 
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
        frq        = self.freq
        lamb       = self.lamb
        phis       = self.phis
        qE0LT      = self.qE0LT
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)

        max_r  = 0.05              # max radial excursion [m]
        r      = sqrt(x**2+y**2)   # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(max_r,ID='_PYO_G:base_map')
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = I0(Kr)            # bessel function I0
        i1     = I1(Kr)            # bessel function I1

        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)

        # ref Teilchen
        wRo = wRi + self.deltaW                           # ref Teilchen energy (O)
 
        # Teilchen
        converter   = WConverter(wRi,frq)
        phiin       = converter.zToDphi(z) + phis         # Teilchen phase (I)
        wo_wi       = qE0LT*i0*cos(phiin)                 # energy kick (Shislo 4.2.3)
        wi          =  converter.Dp2pToW(zp) + wRi        # Teilchen energy (I) dp/p --> dT
        wo          = wi + wo_wi                          # Teilchen energy (O)   
        dw          = wo - wRo                            # Differenz der energy kicks von Teilchen und ref Teilchen (entspricht delta**2)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        """ !!!IMPORTANT!!! SOLL particle (O) """        
        particleRo = copy(particleRi)(tkin = wRo)
        betao      = particleRo.beta
        gammao     = particleRo.gamma
        gbo        = particleRo.gamma_beta

        zo         = betao/betai*z                     # z (O) (4.2.5) A.Shishlo/J.Holmes
        zpo        = converter.DWToDp2p(dw)            # dW --> dp/p (O)

        # if 0: print('z ',z,'zpf ',zpf)

        factor = qE0LT/(m0c2*gbi*gbo)*i1               # common factor
        if r > 0.:
            xp  = gbi/gbo*xp - x/r*factor*sin(phiin)   # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbo*yp - y/r*factor*sin(phiin)
        elif r == 0.:
            xp  = gbi/gbo*xp
            yp  = gbi/gbo*yp

        T = wRo
        f_track = NP.array([x, xp, y, yp, zo, zpo, T, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        # if 0:
        #     arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        #     arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')

        """ the parent reads these attributes below """
        self._particlef = particleRo

        return f_track
class RFC(I):
    """ 
    Rf cavity as product D*Kick*D (DKD-model)
    """
    def __init__(self,
                EzAvg    = 1.,
                label    = 'RFC',
                PhiSoll  = radians(-30.),
                fRF      = 800.e6,
                gap      = 0.024,
                aperture = PARAMS['aperture'],
                dWf      = FLAGS['dWf'],
                length   = 0.,
                mapping  = 't3d',
                SFdata   = None,
                particle = PARAMS['sollteilchen'],
                position = [0, 0, 0],
                next     = None,
                prev     = None):
        if length > 0:
            self.length = length
        else:
            self.length = gap   # die ideale pillbox
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self._EzAvg   = EzAvg*dWf
        self.phis     = PhiSoll
        self.freq     = fRF
        self.gap      = gap
        self.E0L      = self._EzAvg*self.gap
        self.dWf      = dWf
        self.mapping  = mapping
        self.SFdata   = SFdata 
        self._ttf     = None
        self._deltaW  = None
        
        self['viseo']   = 0.25
        
        if self.mapping != 'dyn':
            # ---> DKD models <--- uses RFG for t3d, simple, base, ttf mappings
            dri   = D(length=0.5*self.length, particle=self.particle, aperture=self.aperture)
            drf   = D(length=0.5*self.length, particle=self.particle, aperture=self.aperture)
            kick  = RFG(EzAvg     = self._EzAvg,
                        PhiSoll   = self.phis,
                        fRF       = self.freq,
                        gap       = self.gap,
                        aperture  = self.aperture,
                        dWf       = self.dWf,
                        mapping   = self.mapping,
                        SFdata    = self.SFdata,
                        particle  = self.particle)
            self.triplet = (dri, kick, drf)
            self._deltaW = kick.deltaW
            tkin_f       = self.particle.tkin + self._deltaW   # tkin after acc. gap
            # UPDATE energy for downstream drift after gap
            drf.adjust_energy(tkin_f)
            self._ttf = self._deltaW/(self.E0L*cos(self.phis)) if self.dWf == 1 else 1.
            # in case off ... (really needed ?)
            self.matrix = NP.dot(drf.matrix,NP.dot(kick.matrix,dri.matrix))
            self._particlef = kick.particlef
            DEBUG_OFF("det[RFC.matrix] = {}".format((NP.linalg.det(self.matrix))))
        elif self.mapping == 'dyn':
            # DYNAC gap model with SF-data (E.Tanke, S.Valero)
            # This is no DKD-model. 
            cav             = _DYN_G(self)
            self.triplet    = (cav,)
            self._deltaW    = cav.deltaW
            self._particlef = cav.particlef
            self._ttf       = cav.ttf
            # ACHTUNG! _DYN_G has no matrix, so use DKD with _T3D_G-matrix instead
            dri   = D(length=0.5*self.length, particle=self.particle, aperture=self.aperture)
            drf   = D(length=0.5*self.length, particle=self.particle, aperture=self.aperture)
            kick  = _T3D_G(self)
            self.matrix = NP.dot(drf.matrix,NP.dot(kick.matrix,dri.matrix))
            # correct energy increase in Node.matrix
            self.matrix[Ktp.T,Ktp.dT] = self._deltaW

            
    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(
                    EzAvg         = self.EzAvg,
                    label         = self.label,
                    PhiSoll       = self.phis,
                    fRF           = self.freq,
                    gap           = self.gap,
                    aperture      = self.aperture,
                    dWf           = self.dWf,
                    length        = self.length,
                    mapping       = self.mapping,
                    SFdata        = self.SFdata,
                    particle      = self.particle(tkin),
                    position      = self.position,
                    next          = self.next,
                    prev          = self.prev)
        self._params = _params
        return self

    def map(self,i_track):
        track = copy(i_track)
        for node in iter(self.triplet):
            f_track = node.map(track)
            track = f_track
        DEBUG_OFF('rfc-map {}'.format(f_track))
        return f_track

    def soll_map(self,i_track):
        si,sm,sf = self.position
        track = copy(i_track)
        for node in iter(self.triplet):
            f_track = node.soll_map(track)
            track = f_track
        f_track[SKOO] += sm
        DEBUG_OFF('rfc-soll {}'.format(f_track))
        return f_track

    def make_slices(self, anz = PARAMS['nbslices']):
        slices = []
        if len(self.triplet) == 3:
            dri,kick,drf = self.triplet
            anz = anz/2 if anz != 0 else 1
            dri_slices = dri.make_slices(anz=anz)
            drf_slices = drf.make_slices(anz=anz)
            slices += dri_slices
            slices.append(kick)
            slices += drf_slices
        else:
            slices.append(self)
        return slices

    @property
    def ttf(self):
        return self._ttf
    @property
    def lamb(self):
        return PARAMS['clight']/self.freq
    @property
    def deltaW(self):
        return self._deltaW
    @property
    def EzPeak(self):
        return self['EzPeak']
    @property
    def EzAvg(self):
        return self._EzAvg
    @EzAvg.setter
    def EzAvg(self,value):
        self._EzAvg = self['EzAvg'] = value
    @property
    def particlef(self):
        return self._particlef
class SIXD(D):
    """ 
    Drift with Sixtrack mapping (experimental!) 
    """
    def __init__(self, label="Dsix", particle=PARAMS['sollteilchen'], position=[0., 0., 0.], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self['viseo'] = 0.
        self.off_soll = copy(self.particle) # !!!IMPORTANT!!!

    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(
            label      = self.label, 
            particle   = self.particle(tkin), 
            position   = self.position, 
            length     = self.length,
            aperture   = self.aperture,
            next       = self.next,
            prev       = self.prev)
        self._params   = _params
        return self

    def shorten(self, length):
        shortSIXD = SIXD(label=self.label, particle=self.particle, position=self.position, length=length, aperture=self.aperture)
        shortSIXD._params = self._params
        return shortSIXD

    def map(self, i_track):
        def fpsigma(psigma, soll):
            beta0    = soll.beta
            E0       = soll.e
            m0c2     = soll.e0
            res      = (1+beta0**2*psigma)**2-(m0c2/E0)**2
            res      = sqrt(res)/beta0-1.
            return res

        def einsplusfpsigma(psigma, soll):
            return 1.+fpsigma(psigma, soll)

        def t3d2six(i_track):
            """ Conversion Trace3D ==> Ripken-Schmidt (sixtrack) """
            soll     = self.particle
            x        = i_track[XKOO]       # [0]
            xp       = i_track[XPKOO]      # [1]
            y        = i_track[YKOO]       # [2]
            yp       = i_track[YPKOO]      # [3]
            z        = i_track[ZKOO]       # [4] z
            dp2p     = i_track[ZPKOO]      # [5] dp/p
            T        = i_track[EKOO]       # [6] kinetic energy SOLL
            s        = i_track[SKOO]       # [8] position SOLL

            E0       = soll.e
            beta0    = soll.beta
            p0       = soll.p             # cp SOLL [MeV]
            m0c2     = soll.e0
            p        = p0/(1.-dp2p)
            E        = sqrt(p**2+m0c2**2) # E aus dp2p und p0
            tkin     = E-m0c2
            particle = self.off_soll(tkin = tkin)
            gb       = particle.gamma_beta
            beta     = particle.beta

            px       = gb*m0c2/E0*xp
            py       = gb*m0c2/E0*yp
            sigma    = z
            try:
                psigma = ((beta0/beta/(1.-dp2p))-1.)/beta0**2
            except Exception as ex:
                print('(dp2p, beta, beta0)', (dp2p, beta, beta0))
                print('in t3d2six(): bad psigma')
                sys.exit(1)
            f_track  = NP.array([x, px, y, py, sigma, psigma, T, 1., s, 1.])
            return f_track

        def six2t3d(i_track):
            """ Conversion Ripken-Schmidt (sixtrack) ==> Trace3D """
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
            particle = self.off_soll(tkin = tkin)
            beta     = particle.beta
            gb       = particle.gamma_beta

            xp       = px/(gb*m0c2/E0)
            yp       = py/(gb*m0c2/E0)
            z        = sigma
            dp2p     = 1.-beta0/beta/(1.+beta0**2*psigma)
            f_track  = NP.array([x, xp, y, yp, z, dp2p, T, 1., s, 1.])
            return f_track

        def rps_map(i_track, l):
            """ Ripken-Schmidt sixtrack map """
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
            particle = self.off_soll(tkin = tkin)
            beta     = particle.beta

            xf       = xi + pxi/einsplusfpsigma(psigmai, soll)*l
            pxf      = pxi
            yf       = yi + pyi/einsplusfpsigma(psigmai, soll)*l
            pyf      = pyi
            sigmaf   = sigmai + (1.-(beta0/beta)*(1.+0.5*(pxi**2+pyi**2)/einsplusfpsigma(psigmai, soll)**2))*l
            psigmaf  = psigmai
            f_track  = NP.array([xf, pxf, yf, pyf, sigmaf, psigmaf, T, 1., s, 1.])
            return f_track

        # body map
        f_track     = t3d2six(i_track)
        DEBUG_OFF('t3d-->six\n{}'.format(f_track))
        f_track     = rps_map(f_track, self.length)
        DEBUG_OFF('SIXD.map\n{}'.format(f_track))
        f_track     = six2t3d(f_track)
        DEBUG_OFF('six-->t3d\n{}'.format(f_track))
        # nicht vergessen! adjust total lattice length
        f_track[SKOO] += self.length
        return f_track
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
        p = Particle(50.)
        Dnode = D("Drift",length=l,aperture=5.)
        self.assertEqual(Dnode.length,10.)
        g = p.gamma
        d_matrix = NP.eye(MDIM,MDIM)
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,LKOO] = l
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
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,LKOO] = l
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"shorten")
        self.assertTrue(Dnode.label == "Drift")
        self.assertTrue(Dnode.type == "D")
        self.assertEqual(Dnode.particle.tkin,p.tkin,"tkin")
    def test_Particle_and_K(self):
        print("\b----------------------------------------test_Particle_and_K")
        p = PARAMS['sollteilchen']
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
            k0w  = sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[x,x]  = cos(dphi);        mx[x,xp]  = sin(dphi)/k0w       
                mx[xp,x] = -k0w*sin(dphi);   mx[xp,xp] = mx[x,x]
                mx[y,y]  = cosh(dphi);       mx[y,yp]  = sinh(dphi)/k0w        
                mx[yp,y] = +k0w*sinh(dphi);  mx[yp,yp] = mx[y,y]
            else:
                l = length
                f = 1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = Proton()
        QFnode = QF("QFoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
        self.assertEqual(QFnode.viseo,+0.5,"viseo")
        self.assertTrue(QFnode.matrix[XPKOO,XKOO] < 0.)
        self.assertTrue(QFnode.matrix[YPKOO,YKOO] > 0.)

        p = Proton(100.)
        QFnode = QFnode.adjust_energy(100.)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.adjust_energy")

        length = 0.5
        QFnode = QFnode.shorten(length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.shorten")

        gradient = 1.0; length = 0.15; p = Proton(100.)
        QFnode = QF("QFfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
    def test_QD_Node(self): 
        print("\b----------------------------------------test_QD_Node")
        def matrix(gradient,length,particle,thin):
            x=0; xp=1; y=2; yp=3; z=4; zp=5; T=6; dT=7; S=8; dS=9
            k02  = K(gradient,particle)
            k0w  = sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[y,y]  = cos(dphi);        mx[y,yp]  = sin(dphi)/k0w        
                mx[yp,y] = -k0w*sin(dphi);   mx[yp,yp] = mx[y,y]
                mx[x,x]  = cosh(dphi);       mx[x,xp]  = sinh(dphi)/k0w        
                mx[xp,x] = +k0w*sinh(dphi);  mx[xp,xp] = mx[x,x]
            else:
                l = length
                f = -1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = Proton()
        QDnode = QD("QDfoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
        self.assertTrue(QDnode.matrix[XPKOO,XKOO] > 0.)
        self.assertTrue(QDnode.matrix[YPKOO,YKOO] < 0.)

        length = 0.5
        QDnode = QDnode.shorten(length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"QD.shorten")

        self.assertEqual(QDnode['label'],'QDfoc',"__getitem__")
        self.assertEqual(QDnode['type'],'QD',"__getitem__")
        self.assertEqual(QDnode.viseo,-0.5,"viseo")
        self.assertEqual(QDnode['viseo'],-0.5,"__getitem__")
        QDnode['viseo'] = 0.75
        self.assertEqual(QDnode['viseo'],0.75,"__setitem__")

        gradient = 1.0; length = 0.15; p = Proton(100.)
        QDnode = QD("QDfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
    def test_Wille_Nodes(self):
        print("\b----------------------------------------test_Wille_Nodes")
        kqf  = -1.20    # [1/m**2]
        kqd  = -kqf
        lqf  = 0.20      # [m]
        lqd  = 0.40      # [m]
        rhob = 3.8197   # [m]
        lb   = 1.50       #[m]
        phib = 11.25    # [deg]
        ld   = 0.55       # [m]
        p    = Proton(tkin=100.)
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
        mqf[SKOO,SKOO]  = 1.;         mqf[SKOO,LKOO]   = lqf     # s
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
        mqd[SKOO,SKOO]  = 1.;         mqd[SKOO,LKOO]   = lqd     # s
        # print("\nmqd");print(mqd);print("MQD.matrix");print(MQD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mqd[i,j],MQD.matrix[i,j],msg="MQD",delta=1e-4)
        """ Wille's Dipole = MB """
        MB = SD("MB",2*radians(phib),rhob,particle=p)
        mb = NP.eye(MDIM,MDIM)
        mb[XKOO,XKOO]  = 0.9239;      mb[XKOO,XPKOO]   = 1.4617;        mb[XKOO,ZPKOO] = 0.2908
        mb[XPKOO,XKOO] = -0.1002;     mb[XPKOO,XPKOO]  = mb[XKOO,XKOO]; mb[XPKOO,ZPKOO] = 0.3827
        mb[YKOO,YKOO]  = 1.0;         mb[YKOO,YPKOO]   = 1.5
        mb[YPKOO,YKOO] = -0.0;        mb[YPKOO,YPKOO]  = mb[YKOO,YKOO]
        mb[ZKOO,XKOO]  = -0.3827;     mb[ZKOO,XPKOO]   = -0.2908;       mb[ZKOO,ZPKOO] = 1.1867
        mb[EKOO,EKOO]  = 1.;          mb[EKOO,DEKOO]   = 0.     # tkin
        mb[SKOO,SKOO]  = 1.;          mb[SKOO,LKOO]    = lb     # s
        # print("\nmb");print(mb);print("MB.matrix");print(MB.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mb[i,j],MB.matrix[i,j],msg="MB",delta=1e-4)
        """ Wille's Kantenfokusierung = MEB """
        MEB = _wedge(radians(phib),rhob, t3d_wedge=False)
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
        md[SKOO,LKOO] = ld
        # print("\nmd");print(md);print("MD.matrix");print(MD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(md[i,j],MD.matrix[i,j],msg="MD",delta=1e-4)
        """ Wille's gesamte Zelle mit Sektordipol und Kantenfokussierung = MZ """
        MZ = MQF*MD*MEB*MB*MEB*MD*MQD*MD*MEB*MB*MEB*MD*MQF
        mz = NP.eye(MDIM,MDIM)
        mz[XKOO,XKOO]  = 0.0808;   mz[XKOO,XPKOO]  = 9.7855;        mz[XKOO,ZPKOO]  = 3.1424
        mz[XPKOO,XKOO] = -0.1015;  mz[XPKOO,XPKOO] = mz[XKOO,XKOO]; mz[XPKOO,ZPKOO] = 0.3471
        mz[YKOO,YKOO]  = -0.4114;  mz[YKOO,YPKOO]  = 1.1280
        mz[YPKOO,YKOO] = -0.7365;  mz[YPKOO,YPKOO] = mz[YKOO,YKOO]
        mz[ZKOO,XKOO]  = -0.34707; mz[ZKOO,XPKOO]  = -3.1424;       mz[ZKOO,ZPKOO]  = 4.1844
        mz[SKOO,LKOO]  = 6.0
        # print("\nmz");print(mz);print("MZ.matrix");print(MZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MZ.matrix[i,j],msg="MZ",delta=1e-4)
        """ Wille's gesamte Zelle mit RD-class (same as MZ?) """
        WEDGE = MEB
        MRD   = RD("RD",2*radians(phib),rhob,WEDGE,particle=p)
        MRDZ  = MQF*MD*MRD*MD*MQD*MD*MRD*MD*MQF
        # print("\nmz");print(mz);print("MRDZ.matrix");print(MRDZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MRDZ.matrix[i,j],msg="MRDZ",delta=1e-4)
    def test_SD_and_RD_Nodes(self):
        print("\b----------------------------------------test_SD_and_RD_Nodes")
        rhob   = 3.8197   # [m]
        phib   = 11.25    # [deg]
        p      = Proton(tkin=100.)
        """ slice and recombine SD """
        sd     = SD("SB",2*radians(phib),rhob,particle=p)
        slices = sd.make_slices(anz=3)
        # for item in slices: print(item.matrix); print()
        sdx = slices[0]                                 # combine the slices again
        for i in range(1,len(slices)):
            sdx = sdx * slices[i]
        # print("\nsd"); print(sd.matrix); print("\nsdx"); print(sdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(sd.matrix[i,j],sdx.matrix[i,j])
        """ slice and recombine RD """
        wedge =  _wedge(radians(phib),rhob,t3d_wedge=True)
        rd = RD("RD",2*radians(phib),rhob,wedge,particle=p)
        slices = rd.make_slices(anz=3)
        # for item in slices: print(item.matrix); print()
        rdx = slices[0]
        for i in range(1,len(slices)):
            rdx = rdx * slices[i]
        # print("\nrd"); print(sd.matrix); print("\nrdx"); print(sdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(rd.matrix[i,j],rdx.matrix[i,j])
    def test_MRK_Node(self):
        print("\b----------------------------------------test_MRK_Node")
        class Agent(object):
            counter = 1
            def __init__(self):
                self.counter = Agent.counter
                Agent.counter += 1
            def do_action(self):
                print(F"Agent # {self.counter} here!")
        a = Agent()
        mrk = MRK("MRK",[a,a])
        mrk.add(Agent())
        mrk.add(Agent())
        mrk.add(Agent())
        # mrk.do_actions()
        self.assertEqual(len(mrk.agents),5)
        mrk.adjust_energy(75.)
        self.assertEqual(len(mrk.agents),5)
    def test_GAP_Node(self):
        print("\b----------------------------------------test_GAP_Node")
        EzAvg   = 2.1             #[MV/m]
        phisoll = radians(-30.)
        gap     = 0.022           #[m]
        freq    = 816.e6          #[Hz]
        gap = GAP("GAP",EzAvg,phisoll,gap,freq,dWf=1) # tkin=500 default
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.03766,delta=1.e-4)
        gap = gap.adjust_energy(6.) # tkin=6
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.023817,delta=1.e-4)
    def test_RFG_Node(self):
        print("\b----------------------------------------test_RFG_Node with _T3D_G map")
        EzAvg = 2.1; phisoll = radians(-30.); gap = 0.044; freq = 816.e6
        rfg = RFG("RFG",EzAvg,phisoll,gap,freq)
        # print(rfg.__dict__)
        # print(rfg.gap_model.__dict__)
        # print(rfg.gap_model.matrix)
        self.assertEqual(rfg.mapping,"t3d")
        self.assertEqual(rfg.gap_model.label,"RFG")
        self.assertEqual(rfg.gap_model.type,"_T3D_G")
        self.assertEqual(rfg.gap_model.EzAvg,2.1)
        self.assertEqual(rfg.gap_model.phisoll,radians(-30.))
        self.assertAlmostEqual(rfg.gap_model.deltaW,0.062206,delta=1.e-4)
        self.assertAlmostEqual(rfg.gap_model.matrix[EKOO,DEKOO],0.062206,delta=1.e-4)
        self.assertEqual(rfg.gap_model.length,0.)
        self.assertEqual(rfg.gap_model.matrix[SKOO,LKOO],0.)
        rfg = rfg.adjust_energy(250.)
        self.assertAlmostEqual(rfg.gap_model.deltaW,0.0751,delta=1.e-4)
        self.assertAlmostEqual(rfg.gap_model.matrix[EKOO,DEKOO],0.0751,delta=1.e-4)

if __name__ == '__main__':
    FLAGS['verbose'] = 3
    unittest.main()