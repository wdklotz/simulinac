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
# TODO adjust_energy() and shorten() return new objects. make sure object properties like links are correctly passed.
import sys
from math import sqrt, sinh, cosh, sin, cos, tan, modf, pi, radians, degrees, ceil
from copy import copy, deepcopy
import numpy as NP
import pprint, inspect
import unittest

from setutil import PARAMS, FLAGS, Particle
from setutil import WConverter, dictprnt, objprnt, Proton, Electron
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO, MDIM
from setutil import dBdxprot, scalek0prot, k0prot, I0, I1, arrprnt, Ktp
from Ez0     import SFdata
# from DynacG  import _DYN_G

def PRINT_PRETTY(obj=None):
    file = inspect.stack()[0].filename
    print(F'DEBUG_ON[{file}] ==> ',end="")
    if obj != None: pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj=None):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

# used about everywhere
twopi = 2.*pi
# numpy pretty printing
NP.set_printoptions(linewidth = 132, formatter = {'float': '{:>8.5g}'.format})

class OutOfRadialBoundEx(Exception):
    def __init__(self,max_r,ID=''):
        self.message = "OutOfRadialBoundEx: in '{}' out of {} [cm] max radial excursion.".format(ID,max_r*100.)
""" ------- The mother of all lattice element objects (a.k.a. nodes)# ------ """
class Node(object):
    """ Base class for transfer matrices (linear map)
        ii)  is a dictionary (DictObject base class)
        ii)  each instance holds its copy of the refrence particle (self.particle)
    """
    def __init__(self):
        self.type         = self.__class__.__name__      # self's node type
        # self.map          = self._map_# function ptr to mapping method
        self.particle     = None      # !!!IMPORTANT!!! local copy of the particle object
        self.matrix       = None      # MDIMxMDIM zero matrix used here
        self.position     = None      # [entrance, middle, exit]
        self.length       = None
        self.label        = None
        self.aperture     = None
        self.next         = None      # link to right Node
        self.prev         = None      # link to jeft Node
        self.viseo        = None      # default - invisible
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
        self.sigxy = sigxy
        return sigmas
class I(Node):
    """  Unity matrix: the unity Node """
    def __init__(self, label='I'):
        super().__init__()
        self.label    = label
        self.matrix   = NP.eye(MDIM,MDIM) # set the NODE's member variable
        self.length   = 0.
class MRK(I):
    """ 
    MaRKer node (a.k.a element): Each marker is parent of an agent that does the specific action.
    The action can be bypassed if the 'maction'-FLAG is False.
    """
    def __init__(self, label, active, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.)):
        super().__init__(label)
        self.active     = active
        self.particle   = copy(particle)
        self.position   = position
        self.viseo      = 4.
    def noaction(self,*args):
        pass
class D(Node):
    """  Trace3D drift space  """
    def __init__(self, label, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), length=0.,aperture=None):
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
        m[SKOO, LKOO]  = self.length # Ds the longitudinal length increase

    def adjust_energy(self, tkin):
        adjusted = D(self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture)
        return adjusted
    def shorten(self, length):
        shortend =  D(self.label, particle=self.particle, position=(0.,0.,0.), length=length, aperture=self.aperture)
        return shortend
class DKD(D):
    """  Trace3D drift spacer for Drift_Kick-Drift sandwich (for use with DYNAC CAVNUM cavities) """
    def __init__(self, label, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), length=0.,aperture=None):
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
    def __init__(self, label, grad, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), length=0., aperture=None):
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
        k0w   = sqrt(self.k02)
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
    def __init__(self, label, grad, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), length=0., aperture=None):
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
    def __init__(self, label, alpha, rho, n=0, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
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
        alpha = radians(self.alpha)              # [rad]
        rho   = self.rho                 # [m]
        Ds     = abs(rho*alpha)     # [m]
        h = alpha/Ds      # [1/m]
        kx = sqrt((1-self.n)*h**2)       # [1/m]
        ky = sqrt(self.n*h**2)
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
        m[SKOO, LKOO]  = self.length  # length increase
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
    def __init__(self, label, alpha, rho, wedge, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
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
        beta = radians(self.kwinkel)    
        g = 0.050    # fixed gap of 50 mm assumed
        K1=0.45
        K2=2.8
        psi  = K1*g/rho*((1+sin(beta)**2)/cos(beta))*(1-K1*K2*(g/rho)*tan(beta)) if self.t3d_wedge else 0.
        mxpx = tan(beta)/rho
        mypy = tan(beta-psi)/rho
        mx   = NP.eye(MDIM,MDIM)
        mx[XPKOO, XKOO] = mxpx
        mx[YPKOO, YKOO] = -mypy
        self.matrix = mx
class GAP(Node):
    """ Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    ... nicht sehr nuetzlich: produziert keine long. Dynamik wie Trace3D RFG!  """
    def __init__(self, label, EzAvg, phisoll, gap, freq, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf']):
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
        m[EKOO, DEKOO] = deltaW      # energy increase
        m[SKOO, LKOO]  = self.length # length increase
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
    """  RF-gap of zero length with different kick gap-models """
    def __init__(self, label, EzAvg, phisoll, gap, freq, SFdata=None, particle=Proton(PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=FLAGS['dWf'], mapping='t3d', fieldtab=None):
        super().__init__()
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65) """
            x = gap/(beta*lamb)
            return NP.sinc(x)
        self.particle  = particle
        self.position  = position
        self.length    = 0.
        self.label     = label
        self.aperture  = aperture
        self.viseo     = 0.25
        self.dWf       = dWf                 # dWf=1 wirh acceleration else 0
        self.EzAvg     = EzAvg*self.dWf      # [MV/m] average gap field
        self.phisoll   = phisoll             # [radians] soll phase
        self.freq      = freq                # [Hz]  RF frequenz
        self.omega     = twopi*self.freq
        self.lamb      = PARAMS['clight']/self.freq
        self.gap       = gap                 # [m] rf-gap
        self.mapping   = mapping             # map model
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.E0L       = self.EzAvg*self.gap
        self.qE0LT     = self.E0L*self.ttf
        self.deltaW    = self.E0L*self.ttf*cos(self.phisoll)
        self.particlef = Proton(self.particle.tkin+self.deltaW)
        self.SFdata    = SFdata               # SuperFish data
        self.matrix    = None
        self.map       = None
        self.fieldtab  = fieldtab
        """ dispatching to different gap models """
        if self.mapping   == 't3d' :
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.map = super().map   # delegate mapping to standard matrix multiplication
        elif self.mapping == 'simple':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.map = self.simple_map
        elif self.mapping == 'base':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            self.map = self.base_map_1
        elif self.mapping == 'oxal':
            # self.matrix ->  #OXAL has its own matrix
            self.particlef = None
            # self.map  =>  # OXAL has its own mapping method
        elif self.mapping == 'ttf':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # TTF_G has its own mapping method
        elif self.mapping == 'dyn':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # DYN_G has its own mapping method
        # TODO other mappings not tested TODO
        else:
            print(F"INFO: RFG is a kick-model and does not work with {self.mapping} mapping! Use one of [t3d,simple,base,ttf,oxal].")

    def T3D_matrix(self,ttf, particlei, particlef, E0L, phisoll, lamb, deltaW, length):
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m        = NP.eye(MDIM,MDIM)
        Wav     = particlei.tkin+deltaW/2.   # average tkin
        pav     = Proton(Wav)
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
        m[SKOO, LKOO]  = length
        return m
    def simple_map(self, i_track):
        """ Mapping (i) to (f) in Simplified Matrix Model. (A.Shislo 4.1) """
        xi        = i_track[XKOO]       # [0]
        xpi       = i_track[XPKOO]      # [1]
        yi        = i_track[YKOO]       # [2]
        ypi       = i_track[YPKOO]      # [3]
        zi        = i_track[ZKOO]       # [4] z
        zpi       = i_track[ZPKOO]      # [5] Dp/p
        T         = i_track[EKOO]       # [6] kinetic energy REF
        S         = i_track[SKOO]       # [8] position REF

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
        

        DEBUG_OFF('simple_map: (deltaW,qE0LT,i0,phisoll) = ({},{},{},{})'.format(deltaW,qE0LT,1.,phisoll))

        DW         = deltaW
        WOUT       = WIN + DW
        particlef  = copy(particle)(tkin = WOUT) # !!!IMPORTANT!!!
        betaf      = particlef.beta
        gammaf     = particlef.gamma
        gbf        = particlef.gamma_beta

        # the longitudinal 2x2 map (always linear!) A.Shishlo/J.Holmes (4.1.6-10)
        m11 = gbf/gbi
        m12 = 0.
        m21 = qE0LT*twopi/(lamb*betai)*sin(phisoll)
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
        xpf  = gbi/gbf*xpi - xi * (pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf)) * sin(phisoll) # A.Shishlo/J.Holmes 4.1.11)
        ypf  = gbi/gbf*ypi - yi * (pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf)) * sin(phisoll)

        # track @ out of node
        f_track = NP.array([xi, xpf, yi, ypf, zf, zfp, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        if 0:
            itr = i_track.copy()
            ftr = f_track.copy()
            for i in range(len(f_track)-4):
                itr[i]  = itr[i]*1.e3
                ftr[i]  = ftr[i]*1.e3
            arrprnt(itr, fmt = '{:6.3g},', txt = 'simple_map:i_track:')
            arrprnt(ftr, fmt = '{:6.3g},', txt = 'simple_map:f_track:')
        # TODO 2 lines below still needed?
        # the parent delegates reading these properties from here
        # self._particlef = copy(particle)(particle.tkin + deltaW) # !!!IMPORTANT!!!
        return f_track
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
        r      = sqrt(x**2+y**2)  # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(max_r,ID='_PYO_G:base_map')
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = I0(Kr)                               # bessel function I0
        i1     = I1(Kr)                               # bessel function I1
        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)
        # SOLL
        WIN       = tki                               # energy (i)
        DELTAW    = deltaW                       # energy kick
        WOUT      = WIN + DELTAW                      # energy (f) (4.1.6) A.Shishlo/J.Holmes
        # PARTICLE
        converter = WConverter(WIN,freq)
        # phin      = -z * twopi/(betai*lamb) + phis     # phase (i)  alte methode
        phin      = converter.zToDphi(z) + phisoll          # phase (i)
        deltaW    = qE0LT*i0*cos(phin)                   # energy kick
        # win     = (zp * (gammai+1.)/gammai +1.) * WIN  # energy (i) dp/p --> dT alte methode
        win       =  converter.Dp2pToW(zp) + WIN         # energy (i) dp/p --> dT
        wout      = win + deltaW                         # energy (f)   (4.2.3) A.Shishlo/J.Holmes
        dw        = wout - WOUT                          # d(deltaW)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particlef = Proton(WOUT)       # !!!IMPORTANT!!! SOLL particle (f)
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
            yp  = gbi/gbf*yp - y/r*commonf*sin(phin)  # should be phi-middle
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x, xp, y, yp, z, zpf, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        if 0:
            arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
            arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
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
        deg_phisoll= degrees(phisoll)
        qE0LT      = self.qE0LT
        deltaW     = self.deltaW
        
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
        wRo = wRi + deltaW                           # ref Teilchen energy (O)
 
        # Teilchen
        converter   = WConverter(wRi,freq)
        deg_converter = degrees(converter.zToDphi(z)) 
        phiin       = converter.zToDphi(z) + phisoll 
        deg_phiin = degrees(phiin)        # Teilchen phase (I)
        wo_wi       = qE0LT*i0*cos(phiin)                 # energy kick (Shislo 4.2.3)
        wi          =  converter.Dp2pToW(zp) + wRi        # Teilchen energy (I) dp/p --> dT
        wo          = wi + wo_wi                          # Teilchen energy (O)   
        dw          = wo - wRo                            # Differenz der energy kicks von Teilchen und ref Teilchen (entspricht delta**2)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particleRo = Proton(wRo)
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
        adjusted = RFG(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf,mapping=self.mapping)
        return adjusted

# TODO classes below need unittesting
# TODO D+K+D models not finished for OXAL_C and TTF_C
class RFC(RFG):    #TODO
    """ Rf cavity as product D*Kick*D (DKD-model) """
    def __init__(self, label, EzAvg, phisoll, gap, length, freq, SFdata=None, particle = Proton(PARAMS['injection_energy']), position = (0.,0.,0.), aperture=None, dWf=FLAGS['dWf'], mapping='t3d'):
        super().__init__(label, EzAvg, phisoll, gap, freq, SFdata=None, particle = Proton(PARAMS['injection_energy']), position = (0.,0.,0.), aperture=None, dWf=FLAGS['dWf'], mapping='t3d')
        self.length  = length if length >= gap else gap
        self.dri     = DKD(label="-",particle=particle,position=position,length=self.length/2.,aperture=aperture)
        self.drf     = DKD(label="-",particle=particle,position=position,length=self.length/2.,aperture=aperture)
        self.triplet = (dri,self,drf)
        # UPDATE energy for downstream drift after gap
        tkin_f = self.particle.tkin + self._deltaW   # tkin after acc. gap
        drf.adjust_energy(tkin_f)
        self.matrix = NP.dot(drf.matrix,NP.dot(self.matrix,dri.matrix))
        # DEBUG_OFF("det[RFC.matrix] = {}".format((NP.linalg.det(self.matrix))))
    def adjust_energy(self, tkin):
        adjusted = RFC(self.label,self.EzAvg,self.phisoll,self.gap,self.length,self.freq,SFdata=self.SFdata,particle=Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf,mapping=self.mapping)
        return adjusted
    def map(self,i_track):
        track = copy(i_track)
        for node in iter(self.triplet):
            f_track = node.map(track)
            track = f_track
        DEBUG_OFF('rfc-map {}'.format(f_track))
        return f_track
    def make_slices(self, anz=0):
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
        gradient = 3.; length = 1.; p = Proton(80.)
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
        p = Proton(50.)
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
        p = Proton(50.)
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
        gradient = 3.; length = 1.; p = Proton(80.)
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
        gradient = 3.; length = 1.; p = Proton(80.)
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

        gradient = 1.0; length = 0.15; p = Proton(100.)
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
        p    = Proton(100.)
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
        MB = SD("MB",2*phib,rhob,particle=p)
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
        md[SKOO,LKOO] = ld
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
        mz[SKOO,LKOO]  = 6.0
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
        p      = Proton(100.)

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
        p100      = Proton(100.)
        p50       = Proton(50.)

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
    def test_MRK_Node(self):
        print("\b----------------------------------------test_MRK_Node")
        FLAGS['maction'] = True          # call marker actions
        class Agent(object):
            counter = 1
            def __init__(self):
                self.counter = Agent.counter
                Agent.counter += 1
            def do_action(self):
                return(F"Agent # {self.counter} here!")
        mrks = [MRK(f'marker {i}',Agent(),True) for i in range(3)]
        for cnt, mrk in enumerate(mrks):
            res = f'{mrk.label} {mrk.actions()}'
            # print(res)
            self.assertEqual(res, f'marker {cnt} Agent # {cnt+1} here!')
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
    def test_RFG_Node_with_T3D_G_map(self):
        print("\b----------------------------------------test_RFG_Node_with_T3D_G_map")
        EzAvg = 2.1; phisoll = radians(-30.); gap = 0.044; freq = 816.e6
        rfg = RFG("RFG",EzAvg,phisoll,gap,freq)
        self.assertEqual(rfg.mapping,"t3d")
        self.assertEqual(rfg.label,"RFG")
        self.assertEqual(rfg.EzAvg,2.1)
        self.assertEqual(rfg.phisoll,radians(-30.))
        self.assertAlmostEqual(rfg.deltaW,0.062206,delta=1.e-4)
        self.assertAlmostEqual(rfg.matrix[EKOO,DEKOO],0.062206,delta=1.e-4)
        self.assertEqual(rfg.length,0.)
        self.assertEqual(rfg.matrix[SKOO,LKOO],0.)
        rfg = rfg.adjust_energy(250.)
        self.assertAlmostEqual(rfg.deltaW,0.0751,delta=1.e-4)
        self.assertAlmostEqual(rfg.matrix[EKOO,DEKOO],0.0751,delta=1.e-4)

if __name__ == '__main__':
    FLAGS['verbose'] = 3
    unittest.main()