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

def PRINT_PRETTY(obj):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')
DEBUG_MAP    = False
DEBUG_PYO_G  = False

from setutil import wille, PARAMS, FLAGS
from setutil import WConverter, dictprnt, objprnt, Proton, Electron
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO, MDIM
from setutil import dBdxprot, scalek0prot, k0prot, I0, I1, arrprnt, Ktp
from Ez0     import SFdata
from TTFG    import _TTF_G
from DynacG  import _DYN_G
from OXAL    import _OXAL

twopi = 2.*pi     # used about everywhere

# numpy pretty printing
NP.set_printoptions(linewidth = 132, formatter = {'float': '{:>8.5g}'.format})

#------- The mother of all lattice elements (a.k.a. matrices)
class _Node(object):
    """ Base class for transfer matrices (linear map)
        ii)  is a dictionary (DictObject base class)
        ii)  each instance holds its copy of the refrence particle (self.particle)
    """
    def __init__(self, label='', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=None, next=None, prev=None):
        # DictObject.__init__(self)
        self.matrix    = NP.zeros([MDIM,MDIM])   # MDIMxMDIM zero matrix used here
        # !!!IMPORTANT!!! local copy of the particle object
        self.particle  = copy(particle)
        self.position  = position           # [entrance, middle, exit]
        self.length    = length             # default - thin
        self.label     = label              # default - unlabeled
        self.aperture  = aperture           # default - infinite aperture
        self.next      = next               # right link
        self.prev      = prev               # left link
        self['slice_min'] = 0.001           # default - minimal slice length
        self['viseo']     = 0               # default - invisible
        # a class is a dictionary
        self.type = self.__class__.__name__ # self's node type
        self._params = self.__dict__        # make legacy code compatible

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
        _params = self._params
        self.__init__(label=self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture= self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

    @property
    def twiss(self):
        return self['twiss']

    @property
    def particlef(self):
         # !!!IMPORTANT!!! return a copy with updated energy
        return copy(self.particle)(self.particle.tkin + self.deltaW)

    def __call__(self, n = MDIM, m = MDIM):
        # return upper left n, m submatrix
        return self.matrix[:n, :m]

    def __mul__(self, other):
        """ define the * operator for _Node objects """
        res = _Node()
        if (self.label == ''):
            res.label = other.label
        else:
            res.label = self.label+'*'+other.label
        res.length = self.length + other.length
        """ matrix product """
        res.matrix = NP.dot(self.matrix, other.matrix)
        return res

    def prmatrix(self):
        n  = 42
        nx = 300
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
        raise RuntimeError('_Node:reverse not implemented!')
        sys.exit()

    def inverse(self):
        raise RuntimeError('_Node:inverse not implemented!')
        sys.exit(1)

    def trace(self):
        return self.tracex()+self.tracey

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
        if DEBUG_MAP:
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

class I(_Node):
    """ 
    Unity matrix
    """
    def __init__(self, label='I', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=None, next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self.matrix    = NP.eye(MDIM,MDIM)     # MDIMxMDIM unit matrix used here

class MRK(I):
    """ 
    Marker node (a.k.a element): owns a list of agents that do the actions
    """
    def __init__(self, label='MRK', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=None, next=None, prev=None, agent=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self._agents = []   # the agent-list
        if agent != None: self.add(agent)
        self['viseo'] = 4.
    
    def add(self,agent):
        self._agents += [agent]

    def do_actions(self):
        """ invoke all actions bound to this marker """
        for agent in self._agents:
            agent.do_action()
        
    def adjust_energy(self, tkin):
        _params = self._params
        _agents = self._agents
        self.__init__(label=self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture,  next=self.next, prev=self.prev)
        self._params = _params
        self._agents = _agents
        return self

class D(I):
    """ 
    Trace3D drift space 
    """
    def __init__(self, label='D', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0.,aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        # linear NODE matrix
        m = self.matrix
        g = self.particle.gamma
        m[XKOO, XPKOO] = m[YKOO, YPKOO] = self.length
        m[ZKOO, ZPKOO] = self.length/(g*g)
        m[SKOO, LKOO]  = self.length        # length increase

    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(label=self.label, particle=self.particle(tkin), position=self.position, length=self.length, aperture=self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

    def shorten(self, length):
        shortD =  D(label=self.label, particle=self.particle, length=length, aperture=self.aperture)
        shortD._params = self._params
        return shortD

class QF(D):
    """ 
    Trace3D focussing quad 
    """
    def __init__(self, k0=0., label='QF', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self.k0       = k0         # [m**-2]
        self._mx()
        self['viseo'] = +0.5

    def adjust_energy(self, tkin):
        ki = self.k0
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        kf = ki*cpi/cpf # scale quad strength with new impulse
        _params = self._params
        self.__init__(k0=kf, label=self.label, particle=self.particle, position=self.position, length=self.length, aperture=self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

    def shorten(self, length):
        shortQF = QF(k0=self.k0, label=self.label, particle=self.particle, position=self.position, length=length, aperture=self.aperture)
        shortQF._params = self._params
        return shortQF

    def _mx(self):
        m = self.matrix
        g = self.particle.gamma
        rzz12 = self.length/(g*g)
        kwurz = sqrt(self.k0)
        phi = self.length*kwurz
        # focusing
        cf   = cos(phi)
        sf   = sin(phi)/kwurz
        cfp  = -kwurz*sin(phi)
        sfp  = cf
        # defocusing
        cd  =  cosh(phi)
        sd = sinh(phi)/kwurz
        cdp = kwurz*sinh(phi)
        sdp = cd
        # MDIMxMDIM matrix
        if isinstance(self, QF) and not isinstance(self, QD):
            m[XKOO, XKOO]  = cf; m[XKOO, XPKOO] = sf; m[XPKOO, XKOO] = cfp; m[XPKOO, XPKOO] = sfp
            m[YKOO, YKOO]  = cd; m[YKOO, YPKOO] = sd; m[YPKOO, YKOO] = cdp; m[YPKOO, YPKOO] = sdp
            m[ZKOO, ZPKOO] = rzz12
        elif isinstance(self, QD):
            m[XKOO, XKOO]  = cd; m[XKOO, XPKOO] = sd; m[XPKOO, XKOO] = cdp; m[XPKOO, XPKOO] = sdp
            m[YKOO, YKOO]  = cf; m[YKOO, YPKOO] = sf; m[YPKOO, YKOO] = cfp; m[YPKOO, YPKOO] = sfp
            m[ZKOO, ZPKOO] = rzz12
        else:
            raise RuntimeError('QF: neither QF nor QD! should never happen! - STOP')
            sys.exit(1)
        m[SKOO, LKOO]  = self.length # length increase
        return m

class QD(QF):
    """ 
    Trace3D defocussing quad 
    """
    def __init__(self, k0=0., label='QD', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(k0=k0, label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self['viseo'] = -0.5

    def shorten(self, length):
        shortQD = QD(k0=self.k0, length=length, label=self.label, particle=self.particle, position=self.position, aperture=self.aperture)
        shortQD._params = self._params
        return shortQD

class SD(D):
    """ 
    Trace3d sector dipole in x-plane 
    """
    def __init__(self, radius=0., label='SD', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self.radius   = radius
        self._mx()
        self['viseo'] = 0.25

    def adjust_energy(self, tkin):
        ri = self.radius
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        rf = ri*cpf/cpi # scale bending radius with new impulse
        _params = self._params
        self.__init__(radius=rf, label=self.label, particle=self.particle, position=self.position, length=self.length, aperture=self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

    def shorten(self, length):
        shortSD = SD(radius=self.radius, label=self.label, particle=self.particle, position=self.position, length=length, aperture=self.aperture)
        shortSD._params = self._params
        return shortSD

    def _mx(self):
        m = self.matrix
        rho = self.radius
        k = 1./rho
        phi = self.length/rho
        cx = cos(phi) ; sx = sin(phi)
        b = self.particle.beta
        # x,x'-plane
        m[XKOO, XKOO]  = cx;     m[XKOO, XPKOO]   = sx/k;  m[XKOO, ZPKOO]  = rho*(1.-cx)
        m[XPKOO, XKOO] = -sx*k;  m[XPKOO, XPKOO]  = cx;    m[XPKOO, ZPKOO] = sx
        # y,y'-plane
        m[YKOO, YPKOO] = self.length
        # z,z'-plane
        m[ZKOO, XKOO] = -sx;   m[ZKOO, XPKOO] = -rho*(1.-cx);   m[ZKOO, ZPKOO] = rho*sx-self.length*b*b
        m[SKOO, LKOO] = self.length # length increase
        return m

class RD(SD):
    """ 
    Trace3D rectangular dipole x-plane 
    """
    def __init__(self, radius=0., label='RD', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(radius=radius, label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        psi = 0.5*length/radius   # halber Kantenwinkel

        self.wd = _wedge(psi, radius, particle, position)  # wedge
        rd = self.wd * (self * self.wd)
        self.matrix = rd.matrix

    def make_slices(self, anz=PARAMS['nbslices']):
        DEBUG_OFF('RD.make_slices: {} {:8.4f}'.format(self.label, self.length))
        shortSD = self.shorten(self.length/anz)
        slices = [self.wd]          # wedge @ entrance
        for i in range(anz):
            slices.append(shortSD)
        slices.append(self.wd)      # wedge @ exit
        DEBUG_OFF('slices {}'.format(slices))
        return slices

class _wedge(I):
    """ 
    Trace3d dipole wedge x-plane 
    """
    def __init__(self, psi, radius, particle, position):
        super().__init__(label='w', particle=particle, position=position)
        ckp = tan(psi)/radius
        m = self.matrix
        # MDIMxMDIM matrix
        m[XPKOO, XKOO] = +ckp
        m[YPKOO, YKOO] = -ckp

class GAP(I):
    """ 
    Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    ... nicht sehr nuetzlich: produziert keine long. Dynamik wie Trace3D RFG! 
    """
    def __init__(self,
                    EzAvg      = 1.,
                    PhiSoll    = radians(-30.),
                    fRF        = 800.e6,
                    label      = 'GAP',
                    particle   = PARAMS['sollteilchen'],
                    gap        = 0.024,
                    position   = [0, 0, 0],
                    aperture   = None,
                    dWf        = FLAGS['dWf'],
                    next       = None,
                    prev       = None):
        super().__init__(label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev)
        self.EzAvg    = EzAvg*dWf    # [MV/m] av. gap field strength
        self.phis     = PhiSoll      # [radians] soll phase
        self.freq     = fRF          # [Hz]  RF frequenz
        self.gap      = gap          # [m] eff. gap-length
        self.dWf      = dWf

        self['viseo']  = 0.25

        lamb         = PARAMS['clight']/fRF                  # [m] wellenlaenge
        beta         = self.particle.beta                    # beta Einstein
        tr           = self._trtf_(beta, lamb, gap)          # time-transition factor
        E0L          = self.EzAvg*self.gap                   # Spaltspannung
        self.deltaW  = E0L*tr*cos(PhiSoll)                   # delta-W T.Wrangler pp.221
        tkin         = self.particle.tkin
        tk_center    = self.deltaW*0.5 + tkin                # energy in gap center
        part_center  = copy(self.particle)(tk_center)        # !!!IMPORTANT!!! particle @ gap center
        bg           = part_center.gamma_beta                # beta*gamma @ gap center
        matrix       = self.matrix
        m0c2         = self.particle.e0
        # linear matrix
        self._mx_(matrix, m0c2, E0L, tr, PhiSoll, lamb, bg)

    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(
            EzAvg       = self.EzAvg,
            PhiSoll     = self.phis,
            fRF         = self.freq,
            label       = self.label,
            particle    = self.particle(tkin),
            gap         = self.gap,
            position    = self.position,
            dWf         = self.dWf,
            aperture    = self.aperture, 
            next        = self.next, 
            prev        = self.prev)
        self._params = _params
        return self

    def deltaW(self):
        return self.deltaW

    def _trtf_(self, beta, lamb, gap):
        """ tt-factor nach Panofsky (Lapostolle CERN-97-09 pp.65) """
        teta = gap / (beta*lamb)
        res = NP.sinc(teta)/teta
        return res

    def _mx_(self, m, m0c2, E0L, tr, phis, lamb, bg):
        """ cavity nach Dr.Tiede pp.33 """
        cyp = cxp = -pi*E0L*tr*sin(phis)/(m0c2*lamb*bg**3)
        m[XPKOO, XKOO] = cxp
        m[YPKOO, YKOO] = cyp
        m[EKOO, DEKOO] = self.deltaW # energy increase

class RFG(I):
    """ 
    Wrapper to zero length RF kick gap-models
    """
    def __init__(self,
            EzAvg      = 1.,
            label      = 'RFG',
            PhiSoll    = radians(-30.),
            fRF        = 800.e6,
            gap        = 0.024,
            aperture   = None,
            dWf        = FLAGS['dWf'],
            mapping    = 't3d',
            SFdata     = None,  # return of SFdata (SuperFish data)
            particle   = PARAMS['sollteilchen'],
            position   = [0, 0, 0],
            next       = None,
            prev       = None):
        super().__init__(label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev)
        self._EzAvg   = EzAvg*dWf          # [MV/m] average gap field
        self['EzAvg'] = self._EzAvg
        self.phis     = PhiSoll            # [radians] soll phase
        self.freq     = fRF                # [Hz]  RF frequenz
        self.gap      = gap                # [m] rf-gap
        self.dWf      = dWf
        self.mapping  = mapping            # map model
        self.SFdata   = SFdata             # SuperFish data

        self['viseo'] = 0.25
        # makes the T3D matrix default for RFG
        t3d_g          = _T3D_G(self)
        self.gap_model = t3d_g

        """ set switch to gap model """
        if self.mapping == 't3d':
            # Trace3D-matrix and use linear gap-model
            self.gap_model = t3d_g
        elif self.mapping == 'simple' or self.mapping == 'base':
            # PyOrbit gap-models w/o SF-data
            self.gap_model = _PYO_G(self, self.mapping)
        elif self.mapping == 'ttf':
            # 3 point TTF-RF gap-model with SF-data  (A.Shishlo/J.Holmes)
            self.gap_model = _TTF_G(self)
        elif self.mapping == 'oxal':
            # openXAL gap-model with SF-data  (A.Shishlo/J.Holmes)
            self.gap_model = _OXAL(self)
        elif self.mapping == 'dyn':
            # DYNAC gap model with SF-data (E.Tanke, S.Valero)
            # self.gap_model = _DYN_G(self)  not for RFG anymore!
            print("RFG is a kick-model and does not work with 'dyn'-mapping!")
            print("RFG is a kick-model and does not work with 'dyn'-mapping!")
            print("RFG is a kick-model and does not work with 'dyn'-mapping!",flush=True)
        #    exit(1)
            
    def adjust_energy(self, tkin):
        _params = self._params
        self.__init__(
            EzAvg      = self.EzAvg,
            PhiSoll    = self.phis,
            fRF        = self.freq,
            label      = self.label,
            particle   = self.particle(tkin),
            gap        = self.gap,
            position   = self.position,
            aperture   = self.aperture,
            mapping    = self.mapping,
            SFdata     = self.SFdata,
            dWf        = self.dWf,
            next       = self.next,
            prev       = self.prev)
        self._params = _params
        return self

    @property
    def omega(self):
        return twopi*self.freq
    @property
    def lamb(self):
        return PARAMS['clight']/self.freq
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
    def ttf(self):
        """ delegate to gap-model """
        return self.gap_model.ttf
    @property
    def deltaW(self):
        """ delegate to gap-model """
        return self.gap_model.deltaW
    @property
    def particlef(self):
        """ delegate to gap-model """
        #TODO: better don't use particlef - only deltaW
        return self.gap_model.particlef

    def map(self, i_track):
        """ delegate to gap-model """
        return self.gap_model.map(i_track)

    def soll_map(self, i_track):
        """ delegate to gap-model """
        f_track = self.gap_model.soll_map(i_track)
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
            # correct energy increase in _Node.matrix
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

class _PYO_G(object):
    """ 
    PyOrbit zero length RF gap-model (A.Shishlo,Jeff Holmes) 
    """
    def __init__(self, parent, mapping):
        def trtf(lamb, gap, beta):
            """ Transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65) """
            x = gap/(beta*lamb)
            ttf = NP.sinc(x)   # sinc(x) = sin(pi*x)/(pi*x)
            return ttf
        self.particle   = parent.particle
        self.phis       = parent.phis
        self.lamb       = parent.lamb
        self.freq       = parent.freq
        self.mapping    = mapping
        self.E0L        = parent.EzAvg*parent.gap
        self.ttf        = trtf(parent.lamb, parent.gap, parent.particle.beta)
        self.qE0LT      = self.E0L*self.ttf
        # deltaW soll-energy increase Trace3D (same as Shishlo)
        self._deltaW    = self.qE0LT*cos(self.phis)
        self._particlef = copy(self.particle)(self.particle.tkin+self._deltaW)

        # UPDATE linear NODE matrix with deltaW
        parent.matrix[EKOO, DEKOO] = self.deltaW

        if self.mapping == 'simple':
            self.which_map = self.simple_map
        elif self.mapping == 'base':
            self.which_map = self.base_map
            
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

        # the long. 2x2 map (always linear!) A.Shishlo/J.Holmes (4.1.6-10)
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
        if DEBUG_PYO_G:
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

    def base_map(self, i_track):
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
        
        r      = sqrt(x**2+y**2)                      # radial coordinate
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = I0(Kr)                               # bessel function I0
        i1     = I1(Kr)                               # bessel function I1
        # SOLL
        WIN       = tki                               # energy (i)
        DELTAW    = self.deltaW                       # energy kick
        WOUT      = WIN + DELTAW                      # energy (f) (4.1.6) A.Shishlo/J.Holmes
        # PARTICLE
        converter = WConverter(WIN,frq)
    #   phin      = -z * twopi/(betai*lamb) + phis       # phase (i)  alte methode
        phin      = converter.zToDphi(z) + phis          # phase (i)
        deltaW    = qE0LT*i0*cos(phin)                   # energy kick
    #   win       = (zp * (gammai+1.)/gammai +1.) * WIN  # energy (i) dp/p --> dT alte methode
        win       =  converter.Dp2pToW(zp) + WIN         # energy (i) dp/p --> dT
        wout      = win + deltaW                         # energy (f)   (4.2.3) A.Shishlo/J.Holmes
        dw        = wout - WOUT                          # d(deltaW)

        DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particlef = copy(particle)(tkin = WOUT)       # !!!IMPORTANT!!! SOLL particle (f)
        betaf     = particlef.beta
        gammaf    = particlef.gamma
        gbf       = particlef.gamma_beta

        converter = WConverter(WOUT,frq)
        z         = betaf/betai*z                     # z (f) (4.2.5) A.Shishlo/J.Holmes
    #   zpf       = gammaf/(gammaf+1.) * dw/WOUT      # dW --> dp/p (f)  alte methode
        zpf       = converter.DWToDp2p(dw)            # dW --> dp/p (f)

        commonf = qE0LT/(m0c2*gbi*gbf)*i1             # common factor
        if r > 0.:
            xp  = gbi/gbf*xp - x/r*commonf*sin(phin)  # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbf*yp - y/r*commonf*sin(phin)
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x, xp, y, yp, z, zpf, T, 1., S, 1.])

        # for DEBUGGING
        if DEBUG_PYO_G:
            itr = i_track.copy()
            ftr = f_track.copy()
            for i in range(len(f_track)-4):
                itr[i]  = itr[i]*1.e3
                ftr[i]  = ftr[i]*1.e3
            arrprnt(itr, fmt = '{:6.3g},', txt = 'base_map:i_track:')
            arrprnt(ftr, fmt = '{:6.3g},', txt = 'base_map:f_track:')

        # the parent delegates reading these properties from here
        self._particlef = particlef
        return f_track

class _T3D_G(object):
    """ Mapping (i) to (f) in Trace3D zero length RF-Gap Model """
    def __init__(self, parent):
        def trtf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65) """
            x = gap/(beta*lamb)
            ttf = NP.sinc(x)   # sinc(x) = sin(pi*x)/(pi*x)
            return ttf

        def mx(m, ttf, beta, gamma, particlei, particlef, E0L, phis, lamb, deltaW):
            """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
            e0 = particlei.e0
            kz = 2.*pi*E0L*ttf*sin(phis)/(e0*beta*beta*lamb)
            ky = kx = -0.5*kz/(gamma*gamma)
            bgi = particlei.gamma_beta
            bgf = particlef.gamma_beta
            bgi2bgf = bgi/bgf
            m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
            m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
            m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
            # UPDATE linear NODE matrix with deltaW
            m[EKOO, DEKOO] = deltaW
            return

        # from parent
        matrix         = parent.matrix
        particle       = parent.particle
        EzAvg          = parent.EzAvg
        phis           = parent.phis
        gap            = parent.gap
        lamb           = parent.lamb
        position       = parent.position
        # function scope
        beta           = particle.beta
        ttf            = trtf(lamb, gap, beta)       # Panofski
        deltaW         = EzAvg*gap*ttf*cos(phis)     # deltaW energy kick Trace3D
        tkin           = particle.tkin
        tk_center      = deltaW*0.5+tkin             # energy in gap center
        part_center    = copy(particle)(tk_center)   # !!!IMPORTANT!!! particle @ gap center
        beta_mid       = part_center.beta            # beta @ gap cemter
        gamma_mid      = part_center.gamma           # gamma @ gap center
        particlef      = copy(particle)(tkin+deltaW) # !!!IMPORTANT!!! particle @ (f)
        E0L            = EzAvg*gap

        # the linear NODE matrix for rf gaps
        mx(matrix, ttf, beta_mid, gamma_mid, particle, particlef, E0L, phis, lamb, deltaW)

        # class scope
        # parent may delegate reading these properties from here
        self.matrix     = matrix
        self.ttf        = ttf
        self.deltaW     = deltaW
        self.particlef  = particlef
        self.position   = position

    def map(self, i_track):
        """ Mapping from (i) to (f) with linear Trace3D matrix """
        f_track = copy(i_track)
        f_track = NP.dot(self.matrix,f_track)
        DEBUG_OFF('t3d-map {}'.format(f_track))
        return f_track

    def soll_map(self, i_track):
        si,sm,sf = self.position
        f_track = copy(i_track)
        f_track[EKOO] += self.deltaW
        f_track[SKOO] = sm
        DEBUG_OFF('t3d-soll {}'.format(f_track))
        return f_track

class _thin(_Node):
    """ 
    Base class for thin elements implemented as triplet D*Kick*D 
    """
    def __init__(self,label='',  particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=None, next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev)
        self.matrix    = NP.eye(MDIM)     # MDIMxMDIM unit matrix used here

    def make_slices(self, anz = PARAMS['nbslices']):  
        """ 
        Stepping routine through the triplet
        """
        DEBUG_OFF('_thin.make_slices: {} {:8.4f}'.format(self.label, self.length))
        anz1 = int(ceil(anz/2))
        di   = self.triplet[0]
        df   = self.triplet[2]
        kik  = self.triplet[1]
        d1   = di.shorten(di.length/anz1)
        d2   = df.shorten(df.length/anz1)
        slices = []
        for i in range(anz1):
            slices.append(d1)
        slices.append(kik)      # the Kick
        for i in range(anz1):
            slices.append(d2)
        DEBUG_OFF('slices {}'.format(slices))
        return slices

class QFth(_thin):
    """ 
    Thin F-Quad 
    """
    def __init__(self, k0=0., label='QFT', particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev)
        self.k0        = k0
        self['viseo']  = +0.5
        di = D(length  = 0.5*self.length, particle = self.particle)
        df = di
        kick = _kick(self, particle=self.particle, position=self.position, aperture=self.aperture)
        lens = df * (kick * di) # matrix produkt df*kick*di
        self.matrix = lens.matrix
        self.triplet = (di, kick, df)

    def adjust_energy(self, tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin)
        cpf = self.particle.gamma_beta
        ki = self.k0
        kf = ki*cpi/cpf # scale quad strength with new impulse
        _params = self._params
        self.__init__(k0=kf, label=self.label, particle=self.particle, position=self.position, aperture=self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

class _kick(I):
    """ 
    Matrix for thin lens quad 
    """
    def __init__(self, quad, particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=None):
        super().__init__(label='k', particle=particle, position=position, aperture=aperture)
        m = self.matrix
        m[XPKOO, XKOO] = -quad.k0*quad.length
        m[YPKOO, YKOO] = -m[XPKOO, XKOO]

class QDth(QFth):
    """ 
    Thin D-Quad 
    """
    def __init__(self, k0=0., label='QDT', particle=PARAMS['sollteilchen'], position=[0, 0, 0], aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(k0=-k0, label=label, particle=particle, position=position, aperture=aperture, next=next, prev=prev)
        self['viseo']  = -0.5

class QFthx(D):
    """ 
    Thin F-Quad   (express version of QFth) 
    """
    def __init__(self, k0=0., label='QFT', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
        super().__init__(label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self.k0       = k0
        self['viseo'] = +0.5
        L = self.length
        m = self.matrix
        # thin lens quad matrix (by hand calculation)
        m[0, 0]  = 1. - k0*(L**2)/2.
        m[0, 1]  = L - k0*(L**3)/4.
        m[1, 0]  = -k0*L
        m[1, 1]  = m[0, 0]
        m[2, 2]  = 1. + k0*(L**2)/2.
        m[2, 3]  = L + k0*(L**3)/4.
        m[3, 2]  = +k0*L
        m[3, 3]  = m[2, 2]

    def adjust_energy(self, tkin):
        cpi = self.particle.gamma_beta
        self.particle(tkin) # PARTICLE energy adjusted
        cpf = self.particle.gamma_beta
        ki  = self.k0
        kf  = ki*cpi/cpf # scale quad strength with new impulse
        _params = self._params
        self.__init__(k0=kf, label=self.label, particle=self.particle, position=self.position,length=self.length, aperture=self.aperture, next=self.next, prev=self.prev)
        self._params = _params
        return self

    def make_slices(self, anz=PARAMS['nbslices']):
        slices = [self]
        return slices

class QDthx(QFthx):
    """ 
    Thin D-Quad   (express version of QDth) 
    """
    def __init__(self, k0=0., label='QDT', particle=PARAMS['sollteilchen'], position=[0, 0, 0], length=0., aperture=None, next=None, prev=None):
        super().__init__(k0=-k0, label=label, particle=particle, position=position, length=length, aperture=aperture, next=next, prev=prev)
        self['viseo'] = -0.5

class SIXD(D):
    """ 
    Drift with Sixtrack mapping (experimental!) 
    """
    def __init__(self, label="D#", particle=PARAMS['sollteilchen'], position=[0., 0., 0.], length=0., aperture=PARAMS['aperture'], next=None, prev=None):
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

## Utilities 
class Test(_Node):
    def __init__(self, a, b, c, d, e, f, label = 'test'):
        super().__init__()
        self.matrix = NP.array([[a, b, 0., 0., 0., 0., 0., 0., 0., 0.],
                              [c, d, 0., 0., 0., 0., 0., 0., 0., 0.],
                              [0., 0., a, b, 0., 0., 0., 0., 0., 0.],
                              [0., 0., d, e, 0., 0., 0., 0., 0., 0.],
                              [0., 0., 0., 0., a, b, 0., 0., 0., 0.],
                              [0., 0., 0., 0., e, f, 0., 0., 0., 0.],
                              [0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
                              [0., 0., 0., 0., 0., 0., 0., 1., 0., 0.],
                              [0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
                              [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.],
                              ])
        self.label = label

def k0test(gradient = 0., beta = 0., energy = 0.):
    """
        quad strength as function of energy and gradient
        gradient in [Tesla/m]
        energy in [Gev]
        beta [v/c]
    """
    if (gradient !=  0. and energy !=  0. and beta != 0.):
        return 0.2998*gradient/(beta*energy)
    else:
        raise RuntimeError('zero gradient or energy or beta in quad strength!')
        sys.exit(1)
        
def test0():
    print('--------------------------------Test0---')
    print('trivial test 0 ...')
    a = Test(1, 2, 3, 4, 5, 6, label = 'a')
    print(a.string())
    b = Test(1, 1, 1, 1, 1, 1, label = 'b')
    print(b.string())
    print((a*b).string())
    print((b*a).string())
def test1():
    print('--------------------------------Test1---')
    print('trivial test 1 ...')
    i1 = I()
    i2 = i1*i1
    print(i1.string())
    print(i2.string())
def test2():
    print('--------------------------------Test2---')
    print('trivial test 2 ...')
    i1 = I()
    d1 = D(length=10., label='D1')
    print(d1.string())
    print((d1*d1).string())
    print((d1*d1*d1).string())
    print((d1*i1*d1).string())
    d2 = D(length=90, label='D2')
    print(d2.string())
    print((d1*d2).string())
    d3 = D(length=90., label = '')
    print((d2*d3).string())
def test3():
    print('--------------------------------Test3---')
    print('test product of _Node class ...')
    gradient = 1.
    beta     = 0.5
    energy   = 0.2
    print('gradient[Tesla/m] {:.3f}; beta[v/c] {:.3f}; energy[Gev] {:.3f}'.format(gradient, beta, energy))
    k = k0test(gradient = gradient, energy = energy, beta = beta)
    qf = QF(k0 = k, length = 1.1)
    print("QF-->", qf.string())
    # test product of _Node class
    qd = QD(k0 = k, length = 1.2)
    print("QD-->", qd.string())
    print("QF*QD-->", (qf*qd).string())
def test4():
    print('--------------------------------Test4---')
    def doit(elm, anz):
        elm_slices = elm.make_slices(anz = anz)
        print(''.join('{}\n'.format(el) for el in elm_slices))
        elmx = elm_slices[0]
        for slice in elm_slices[1:]:
            elmx = elmx*slice
        print(elmx.string())
        print(elm.string())

    print('test slicing of elements ...')
    gradient = 1.
    beta     = 0.5
    energy   = 0.2
    k = k0test(gradient = gradient, energy = energy, beta = beta)
    # elements
    x = pi
    d      = D(length = x)
    qf     = QF(k0 = k, length = x)
    qd     = QD(k0 = k, length = x)
    sd     = SD(radius = 10., length = x)
    rd     = RD(radius = 10., length = x)
    gp     = GAP()
    rg     = RFG()
    qfth   = QFth(k0 = k)
    qdth   = QDth(k0 = k)
    qfthx  = QFthx(k0 = k, length = x)
    qdthx  = QDthx(k0 = k, length = x)
    rfc    = RFC(length = x)
    sixd   = SIXD(length = x)
    # slicing
    anz = 5
    doit(d, anz)
    doit(sixd, anz)
    doit(qf, anz)
    doit(qd, anz)
    doit(sd, anz)
    doit(rd, anz)
    doit(gp, anz)
    doit(rg, anz)
    doit(qfth, anz)
    doit(qdth, anz)
    doit(qfthx, anz)
    doit(qdthx, anz)
    doit(rfc, anz)
def test5():
    print('--------------------------------Test5---')
    print("K.Wille's Beispiel auf pp.113 Formel (3.200)")
    kqf = wille()['k_quad_f']
    lqf = wille()['length_quad_f']
    kqd = wille()['k_quad_d']
    lqd = wille()['length_quad_d']
    rhob = wille()['bending_radius']
    lb  = wille()['dipole_length']
    ld  = wille()['drift_length']
    # elements
    mqf = QF(k0=kqf, length=lqf, label='QF')
    mqd = QD(k0=kqd, length=lqd, label='QD')
    mb = RD(radius=rhob, length=lb, label='B')
    md = D(length=ld)
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
    kqf =   wille()['k_quad_f']
    lqf =   wille()['length_quad_f']
    kqd =   wille()['k_quad_d']
    lqd =   wille()['length_quad_d']
    rhob =  wille()['bending_radius']
    lb =    wille()['dipole_length']
    ld =    wille()['drift_length']

    # elements
    mqf = QF(k0=kqf, length=lqf, label='QF')
    mqd = QD(k0=kqd, length=lqd, label='QD')
    mb = RD(radius=rhob, length=lb, label='B')
    md = D(length=ld)
    rfc = RFC(length = 4*0.024)

    steps = 13

    # test step_through elements ...
    list = [mqf, mqd, mb, md]
    list = [mqf, rfc]
    for m_anfang in list:
        m_end = I()
        slices = m_anfang.make_slices(anz = steps)
        print('======================================')
        for count, mi in enumerate(slices):
            print('step ', count+1, end = '  ')
            print(mi.string())
            m_end = m_end*mi
        print(m_end, '\n'+m_end.string())
        print(m_anfang, '\n'+m_anfang.string())
def test7():
    print('--------------------------------Test7---')
    print('test Rechteckmagnet ...')
    rhob = wille()['bending_radius']
    lb  = wille()['dipole_length']
    mb  = SD(radius = rhob, length = lb, label = 'SD')
    mr  = RD(radius = rhob, length = lb, label = 'RD')
    print(mb.string())
    print(mr.string())
def test8():
    from lattice import Lattice
    from tracker import track_soll
    print('--------------------------------Test8---')
    print('soll-particle\n'+PARAMS['sollteilchen'].string())
    print('test rf-gap ...')

    gap = GAP()
    lg = Lattice()
    lg.add_element(gap)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'GAP', filter = 'matrix')
    objprnt(gap, 'GAP')
    print(solltrack)
    print('GAP.particle(i)\n'+gap.particle.string())
    print('GAP.particle(f)\n'+gap.particlef.string())

    rfg = RFG()
    lg = Lattice()
    lg.add_element(rfg)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'RFG', filter = 'matrix')
    objprnt(rfg, 'RFG')
    print(solltrack)
    print('RFG.particle(i)\n'+rfg.particle.string())
    print('RFG.particle(f)\n'+rfg.particlef.string())

    rfg = RFG(mapping = 'simple')
    lg = Lattice()
    lg.add_element(rfg)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'RFG', filter = 'matrix')
    objprnt(rfg, 'RFG')
    print(solltrack)
    print('RFG.particle(i)\n'+rfg.particle.string())
    print('RFG.particle(f)\n'+rfg.particlef.string())

    rfg =  RFG(mapping = 'base')
    lg = Lattice()
    lg.add_element(rfg)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'RFG', filter = 'matrix')
    objprnt(rfg, 'RFG')
    print(solltrack)
    print('RFG.particle(i)\n'+rfg.particle.string())
    print('RFG.particle(f)\n'+rfg.particlef.string())

    input_file = 'SF_WDK2g44.TBL'
    EzPeak     = 1.4
    SF_tab    = SFdata(input_file, EzPeak = EzPeak)
    rfg = RFG(mapping = 'ttf', gap = 0.048, SFdata = SF_tab)
    lg = Lattice()
    lg.add_element(rfg)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'RFG', filter = 'matrix')
    objprnt(rfg, 'RFG')
    print(solltrack)
    print('RFG.particle(i)\n'+rfg.particle.string())
    print('RFG.particle(f)\n'+rfg.particlef.string())

    input_file = 'SF_WDK2g44.TBL'
    EzPeak     = 1.4
    SF_tab    = SFdata(input_file, EzPeak = EzPeak)
    rfg = RFG(mapping = 'dyn', gap = 0.048, SFdata = SF_tab)
    lg = Lattice()
    lg.add_element(rfg)
    solltrack = track_soll(lg).as_table()
    # objprnt(rfg, 'RFG', filter = 'matrix')
    objprnt(rfg, 'RFG')
    print(solltrack)
    print('RFG.particle(i)\n'+rfg.particle.string())
    print('RFG.particle(f)\n'+rfg.particlef.string())

    print('test cavity ...')
    cav = RFC()
    lg = Lattice()
    lg.add_element(cav)
    solltrack = track_soll(lg).as_table()
    # objprnt(cav, 'RFC', filter = 'matrix')
    objprnt(cav, 'RFC')
    print(solltrack)
    print('CAV.particle(i)\n'+cav.particle.string())
    print('CAV.particle(f)\n'+cav.particlef.string())
def test9():
    print('--------------------------------Test9---')
    print('test: quad k-faktor and quad scaling ...')
    grad = 45.                              # [T/m] gradient
    tk  = 50.                               # [MeV]  kin. energy
    kq = k0prot(gradient = grad, tkin = tk) # quad strength [1/m**2]
    len = 0.4                               # quad len [m]
    focal = kq*len
    focal = 1./focal  # focal len [m]

    print('soll-particle\n'+PARAMS['sollteilchen'].string())
    print('kq [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(grad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))

    grad = dBdxprot(kq, tk) # quad gradient from k and tkinetic
    print('dB/dz[T/m]\t{:.3f} from dBxprot()'.format(grad))

    mqf = QF(k0 = kq, length = len)
    mqd = QD(k0 = kq, length = len)
    cavity = RFC(
        EzAvg = 1.,
        PhiSoll = radians(-30.),
        fRF = 800.e6)
    print('======================== adjust_energy QF')
    tki = 50.              # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0., 950.]:
        tkf = tki+dt
        k_scaled = scalek0prot(kq, tki, tkf)
        print('(tki, kq) ({},{:.3f}) --> (tkf, k_scaled) ({},{:.3f})'.format(tki, kq, tkf, k_scaled))
        print(mqf.adjust_energy(tkf).string())
        print(mqf.particle.string())
    print('======================== adjust_energy QD')
    tki = 50. #  PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0., 950.]:
        tkf = tki+dt
        k_scaled = scalek0prot(kq, tki, tkf)
        print('(tki, kq) ({},{:.3f}) --> (tkf, k_scaled) ({},{:.3f})'.format(tki, kq, tkf, k_scaled))
        print(mqd.adjust_energy(tkf).string())
        print(mqd.particle.string())
    print('======================== adjust_energy CAV')
    tki = 50.    #  PARAMS['injection_energy']    # [MeV]  kin. energy
    PARAMS['sollteilchen'] = Proton(tki)
    for dt in [0., 950.]:
        tkf = tki+dt
        k_scaled = scalek0prot(kq, tki, tkf)
        print('(tki, kq) ({},{:.3f}) --> (tkf, k_scaled) ({},{:.3f})'.format(tki, kq, tkf, k_scaled))
        print(cavity.adjust_energy(tkf).string())
        print(cavity.particle.string())
def test10():
    print('--------------------------------Test10---')
    print('Particle class test ...')
    dictprnt(PARAMS, text = 'setutil.PARAMS')
    # particle class
    print()
    print(Proton(0.).string())
    print(Proton(50.).string())
    print(Proton(200.).string())
    print(Proton(1.e3).string())
    print()
    print(Electron(0.).string())
    print(Electron(50.).string())
    print(Electron(200.).string())
    print(Electron(1.e3).string())
def test11():
    print('--------------------------------Test11---')
    print('thin lense tests ...')
    k0     =  1.
    length =  2.
    qf     =  QFth(k0 = k0)
    qd     =  QDth(k0 = k0)
    rf     =  RFC(length = length)
    print(qf.string())
    print('soll-particle@QFT\n'+qf.particle.string())
    print(qd.string())
    print('soll-particle@QDT\n'+qf.particle.string())
    print(rf.string())
    print('soll-particle@RFC\n'+qf.particle.string())
    print('---------------- step through ---------------')
    for elm in qf.make_slices(anz = 8):
        print(elm.string())
    print('---------------- step through ---------------')
    for elm in qd.make_slices(anz = 7):
        print(elm.string())
    print('------ RF cavity test & step through --------')
    for elm in rf.make_slices():
        print(elm.string())
def test12():
    print('--------------------------------Test12---')
    print('test12 adjust_energy change ...')
    d = D(length = 99.);              print('id >>', d);     print(d.string())
    d.adjust_energy(tkin = 1000.);    print('id >>', d);     print(d.string())
    qf = QF(k0 = 1.5, length = 0.3);  print('id >>', qf);    print(qf.string())
    qf.adjust_energy(tkin = 200.);    print('id >>', qf);    print(qf.string())
    qd = QD(k0 = 1.5, length = 0.3);  print('id >>', qd);    print(qd.string())
    qd.adjust_energy(tkin = 200.);    print('id >>', qd);    print(qd.string())
    rfc = RFC(length = 1.23);         print('id >>', rfc);   print(rfc.string())
    rfc.adjust_energy(tkin = 200.);   print('id >>', rfc);   print(rfc.string())
def test13():
    print('--------------------------------Test13---')
    print('test SIXD node tracking ...')
    particle = Proton(tkin = 100.)
    l   =  0.05    #[m]
    sixd = SIXD(length = l, particle = particle)
    xi  = yi = 1.e-2
    xpi = ypi = 1.e-2
    z   = 1.e-3
    dp2p = 1.e-2
    i_track = NP.array([xi, xpi, yi, ypi, z, dp2p, 0., 1., 0., 1.])
    f_track = sixd.map(i_track)
    print(i_track)
    print(f_track)
def test14():
    print('--------------------------------Test14---')
    print('test MRO for QF, QD ...')
    print('type(QF.__mro__) =  ', type(QF.__mro__))
    print(''.join('{}\n'.format(el) for el in QF.__mro__))
    print('type(QD.__mro__) =  ', type(QD.__mro__))
    print(''.join('{}\n'.format(el) for el in QD.__mro__))

    qf0 = QF(k0 = 1.0, length = 1.98, label = 'qf0')
    qf1 = QF(k0 = 1.0, length = 1.88, label = 'qf1')
    qf1['viseo'] = 0.6
    qd0=QD(k0 = 1.0, length = 1.78, label = 'qd0')
    qd0['viseo'] = -0.7
    print('qf0.label = ', qf0.label)
    print('qf0["viseo"]= ', qf0['viseo'])
    print('qf1.label = ', qf1.label)
    print('qf1["viseo"]= ', qf1['viseo'])
    print('qd0.label = ', qd0.label)
    print('qd0["viseo"]= ', qd0['viseo'])
    print('\n', qf0.__dict__)
    print('\n', qd0.__dict__)

#------main ----------
if __name__ == '__main__':
    FLAGS['verbose'] = 3
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
    test13()
    test14()
