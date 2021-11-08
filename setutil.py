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
# Python 2 and 3 print compatability
from __future__ import print_function

import sys
import scipy.constants as C
from math import pi,sqrt,sin,cos,radians,degrees,fabs,exp,atan
import logging, pprint
from enum import IntEnum
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import warnings
import time
import pprint

# MDIM: dimension of matrices
MDIM = 10

# new DEBUG facility (replaces old DEBUG_ON,DEBUG_OFF and DEBUG)
def PRINT_PRETTY(obj):
    pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

# Logger
ch        = logging.StreamHandler()     # console handler
ch.setLevel(logging.DEBUG)              # set handler level
formatter = \
    logging.Formatter("%(levelname)s: %(filename)s[%(lineno)d] %(message)s")
ch.setFormatter(formatter)              # set handler's format
logger    = logging.getLogger("logger")
logger.addHandler(ch)                   # add handler to logger

#------  DEFAULT "FLAGS" & "PARAMS" and global dicts
FLAGS  = dict(
        periodic             = False,            # periodic lattice? default
        egf                  = False,            # emittance grow flag default
        sigma                = True,             # beam sizes by sigma-tracking
        KVout                = False,            # print a dictionary of Key-Value pairs, no display
        dWf                  = 1,                # acceleration on/off flag 1=on,0=off
        verbose              = 0,                # print flag default = 0
        express              = False,            # use express version of thin quads
        useaper              = False,            # use aperture check for quads and rf-gaps
        bucket               = False,            # plot bucket
        csTrak               = True,             # plot CS trajectories
        pspace               = False             # plot CS twiss ellipses at entrance
        )
PARAMS = dict(
        clight               = C.c,              # [m/s] const
        elementarladung      = C.e,              # [coulomb] const
        proton_mass          = C.value('proton mass energy equivalent in MeV'),
        electron_mass        = C.value('electron mass energy equivalent in MeV'),
        aperture             = None,             # default aperture = no aperture
        nbwindgs             = 30,               # nboff coil windings
        nbsigma              = 3,                # nboff beam sigmas to stay clear of aperture
        emitx_i              = 2.0e-6,           # [m*rad] Vorgabe emittance entrance
        emity_i              = 2.0e-6,           # [m*rad] Vorgabe emittance entrance
        betax_i              = 2.800,            # [m] Vorgabe twiss beta entrance
        betay_i              = 0.200,            # [m] Vorgabe twiss beta entrance
        alfax_i              = 0.0,              # Vorgabe twiss alpha entrance
        alfay_i              = 0.0,              # Vorgabe twiss alpha entrance
        alfaw_i              = 0.0,              # Vorgabe twiss alpha entrance
        nbslices             = 2,                # default number of slices
        mapset               = frozenset(['t3d','simple','base','ttf','dyn','oxal']), #gap-models
        mapping              = 'base',           # default rf gap-model      
        DT2T                 = 1.e-3,            # default kinetic energy spread  (T a.k.a W)
        warnmx               = 5                 # limit nbof warnings
        )
ELEMENTS = {}
SECTIONS = {}
LATTICE  = []
# (global) SUMMARY: dictionary for summary
SUMMARY = {}

# using enum.IntEnum (since Python 3.4) fuer Koordinatenindizees
# TODO: besser mit namedtupel ?
class Ktp(IntEnum):
    """ Koordinaten fuer track points (1x10)"""
    x  = 0     # x
    xp = 1     # x'
    y  = 2     # y
    yp = 3     # y'
    z  = 4     # z
    zp = 5     # z' = Dp/p
    T  = 6     # T(s) = Ingegral(dT)
    dT = 7     # const 1
    S  = 8     # S = Integral(dS)
    dS = 9     # const 1
    
# for compatability with elder code TODO: replace by namedtupel
XKOO=Ktp.x; XPKOO=Ktp.xp; YKOO=Ktp.y; YPKOO=Ktp.yp; ZKOO=Ktp.z; ZPKOO=Ktp.zp; EKOO=Ktp.T; DEKOO=Ktp.dT; SKOO=Ktp.S; LKOO=Ktp.dS

class Ktw(IntEnum):
    """ Koordinaten fuer twiss vector (1x10) """
    bx = 0      # twiss-beta
    ax = 1      # twiss-alpha
    gx = 2      # twiss-gamma
    by = 3
    ay = 4
    gy = 5
    bz = 6
    az = 7
    gz = 8
    s  = 9      # abszisse for twiss functions

class Twiss(object):
    def __init__(self, beta, alfa, epsi):
        self.epsi  = epsi
        self.beta  = beta
        self. alfa = alfa
        if beta >= 0.:
            self.gamma = (1.+self.alfa**2)/self.beta
        else:
            self.gamma = None
    def __call__(self):
        return (self.beta,self.alfa,self.gamma,self.epsi)   # (beta,alfa,gamma,epsi)
    def y1(self):
        # intersection point(x1,0)
        return (sqrt(self.epsi/self.gamma), 0.)
    def y2(self):
        # maximum point (xmax,x')
        return (sqrt(self.epsi*self.beta),-self.alfa*sqrt(self.epsi/self.beta))
    def y3(self): 
        # maximum point (x,x'max)
        return (-self.alfa*sqrt(self.epsi/self.gamma), sqrt(self.epsi*self.gamma))
    def y4(self):
        # intersection point (0,x4')
        return (0., sqrt(self.epsi/self.beta))
    def sigmaH(self):
        (x,y) = self.y2()
        return x
    def sigmaV(self):
        (x,y) = self.y3()
        return y

class Particle(object):
    """
        A particle class
    """
    def __init__(self,tkin=0.,mass= PARAMS['proton_mass'],name='proton'):
        self._set_self(tkin,mass,name)
    def _set_self(self,tkin,mass,name):
        self.tkin       = tkin                     # kinetic energy [MeV]
        self.e0         = mass                     # rest mass [MeV]
        self.e          = self.e0+self.tkin        # total energy [MeV]
        self.gamma      = self.e/self.e0
        self.lost       = False
        self.track      = None
        try:
            self.beta   = sqrt(1.-1./(self.gamma*self.gamma))
        except ValueError as ex:
            print("Particle's kinetic energy went negative! (tkin[MeV] = {:6.3f})".format(tkin))
            raise ex
            # sys.exit(1)
        self.gamma_beta = self.gamma * self.beta
        self.p          = self.gamma_beta * self.e0   # impulse [Mev]
        self.v          = self.beta * PARAMS['clight']    # velocity [m/s]
        self.brho       = 1.e+6 / PARAMS['clight'] * self.gamma_beta * self.e0 # [T*m]
        self.name       = name
        self.m0c2       = self.e0
        self.m0c3       = self.e0 * PARAMS['clight']
        self.betac      = self.v
        self.E          = self.e
        self.T          = self.tkin
    def string(self):
        headr = ['particle','B*rho[Tm]','Tk[Mev]','p[Mev/c]','gamma','beta','gamma*beta','E[Mev]']
        records = [[
                '{:8s}'.format(self.name),
                '{:8.4f}'.format(self.brho),
                '{:8.4f}'.format(self.tkin),
                '{:8.4f}'.format(self.p),
                '{:8.4f}'.format(self.gamma),
                '{:8.4f}'.format(self.beta),
                '{:8.4f}'.format(self.gamma_beta),
                '{:8.4f}'.format(self.e)
                ]]
        return tblprnt(headr,records)
    def trtf(self,gap,freq):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*freq*gap / (self.beta* PARAMS['clight'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def __call__(self,tkin):  # make callable to change energy
        self._set_self(tkin=tkin,mass=self.e0,name=self.name)
        return self

class Proton(Particle):
    def __init__(self,tkin= 50.):
        super(Proton,self).__init__(tkin=tkin,mass= PARAMS['proton_mass'],name='proton')

class Electron(Particle):
    def __init__(self,tkin= 50.):
        super(Electron,self).__init__(tkin=tkin,mass= PARAMS['electron_mass'],name='electron')

# Sollteichen
PARAMS['sollteilchen'] = Proton()

class WConverter(object):
    """
    Converter to switch between different longitudinal phase space coordinates
    """
    def __init__(self,tk,freq):
        self.pi             = C.pi
        self.lamb           = C.c/freq           # [m]
        self.m0c2           = C.value('proton mass energy equivalent in MeV')
        self.gamma          = 1. + tk/self.m0c2
        self.beta           = sqrt(1.-1./(self.gamma*self.gamma))
        self.b              = self.beta
        self.g              = self.gamma
        self.gb             = self.g*self.b
        self.g2             = self.g**2
        self.bl             = self.b*self.lamb
        self.twopi          = 2.*self.pi

    def DphiToz(self,Dphi):
        """ dela-phi [rad] to z [m] """
        z = -self.bl/self.twopi*Dphi
        return z    # [m]
    def zToDphi(self,z):
        """ z [m] to delta-phi [rad] """
        Dphi = -self.twopi/self.bl*z
        return Dphi    # [rad]

    def wToDp2p(self,w):
        """ w=DW/m0c2(a.k.a. Dgamma) [] to Dp2p [] """
        Dp2p = self.g/(self.g2-1.)*w
        return Dp2p # []
    def Dp2pTow(self,Dp2p):
        """ Dp2p [] to w=DW/m0c2(a.k.a. Dgamma) [] """
        w = (self.g2-1.)/self.g*Dp2p
        return w # []
    def Dp2pToW(self,Dp2p):
        """ Dp2p [] to DW=Dw*m0*c^2 [MeV] """
        W = self.Dp2pTow(Dp2p)*self.m0c2
        return W
    def DWToDp2p(self,DW):
        """ DW [MeV] to Dp2p [] """
        Dp2p = self.wToDp2p(DW/self.m0c2)
        return Dp2p

    def emitwToemitz(self,emitw):
        """ emittance[Dphi(x)w] [rad] to emittance[z(x)Dp2p] [m] """
        emitz = self.lamb/(self.twopi*self.gb)*emitw
        return emitz # [m]
    def emitzToemitw(self,emitz):
        """ emittance[z(x)Dp2p] [m] to emittance[w(x)Dphi] [rad] """
        emitw = (self.twopi*self.gb)/self.lamb*emitz
        return emitw # [rad]

    def betawTobetaz(self,betaw):
        """ beta[w(x)Dphi] [rad] to beta[z(x)Dp2p] [m] """
        betaz = (self.bl)/(self.twopi)*(self.g2-1.)/self.g*betaw
        return betaz # [m]
    def betazTobetaw(self,betaz):
        """ beta[z(x)Dp2p] [m] to beta[w(x)Dphi] [rad] """
        betaw = (self.twopi)/(self.bl)*self.g/(self.g2-1.)*betaz
        return betaw # [rad]
    
    def wtoz(self,args):
        # conversion {Dpi(x)w} to {z(x)Dp2p}
        """ 
        IN: args is tuple = (Delta-phi, w (== dT/m0c2), emittance-w, beta-w)
        OUT:     is tuple = (z, Dp2p, emittance-z, beta=z)
        """
        Dphi, w, emitw, betaw = args
        z = self.DphiToz(Dphi)
        Dp2p = self.wToDp2p(w)
        emitz = self.emitwToemitz(emitw)
        betaz = self.betawTobetaz(betaw)
        return (z,Dp2p,emitz,betaz)

class Functions(object):
    """
    A class to gather function-values (Ordinaten) over a common independent variable (Abszisse)
    """
    def __init__(self,names):
        self._values  = [] # [(abzisse, ordinate-1, ordinate-2, ordinate-3,...)]
        self._points = 0
        nmap = {}
        for i in range(len(names)):
            nmap[names[i]] = i   # {name:i}..... 
        self._nmap = nmap
    def append(self,abszisse,ordinaten):
        value = [abszisse]
        for ord in iter(ordinaten):
            value.append(ord)
        value = tuple(value)
        self._values.append(value)
        self._points += 1
    @property
    def nbpoints(self):
        return self._points
    @property
    def fnames(self):
        return self._nmap
    @property
    def nbfunctions(self):
        return len(self._nmap)-1
    def __getitem__(self,n):
        return self._values[n], self._nmap
    def __call__(self,npoint,nfunc):
        nf = self._nmap[nfunc]
        value = self._values[npoint][nf]
        return value
    pass
    
class SCTainer(object):
    """
    A (singleton) container for objects
    """
    class _singleton_(object):
        def __init__(self):
            self.objects = []
    instance = None
    def __init__(self):
        if not SCTainer.instance:
            SCTainer.instance = self._singleton_()
    def object(self,n):
        return SCTainer.instance.objects[n]
    def objects(self):
        return SCTainer.instance.objects
    def append(self,obj):
        SCTainer.instance.objects.append(obj)
    def len(self):
        return len(SCTainer.instance)

#------------- TimeStamps
class TmStamp(object):
    """
    Utility to place time stamps in the program flow.
    """
    # Remark: works like a sigleton object because it uses only class-attributes and -members.
    tmStamps = []
    tmcnt    = 0
    tmt0     = time.perf_counter()
    @classmethod
    def stamp(self,text):
        self.tmcnt +=1
        t = time.perf_counter()
        dt = t-self.tmt0
    #    self._t0 = t
        s = '{:4}:{:10.5} {:20}'.format(self.tmcnt,dt,text)
        self.tmStamps.append(s)
        return s
    @classmethod
    def as_str(self):
        cntr=0
        str = ''
        for entry in self.tmStamps:
            if cntr%3 == 0:
                str= '{}  {}\n'.format(str,entry)
            else:
                str= '{}  {}'.format(str,entry)
            cntr+=1
        return str
    
## Long. Emittance
def waccept(node):
    """
    Central to calculate longitudinal phase space ellipse parameters nach T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0*Dphi0 = ellipse_area/pi
    IN
        node: the 1st rf-gap at the linac entrance
    """
    if node is not None and FLAGS['dWf'] == 1:
        DT2T      = PARAMS['DT2T']
        E0T       = node.EzAvg*node.ttf  # [MV/m]
        particle  = node.particle
        phis      = node.phis            # [rad]
        lamb      = node.lamb            # [m]
        freq      = node.freq            # [Hz]
        m0c2      = particle.e0          # [MeV]
        gb        = particle.gamma_beta
        beta      = particle.beta
        gamma     = particle.gamma
        tkin      = particle.tkin
        conv      = WConverter(tkin,freq)

        # LARGE amplitude oscillations (T.Wangler pp. 175)
        # w = Dgamma = DW/moc2 is normalized energy deviation
        factor_phis = phis*cos(phis)-sin(phis)
        wmx  = 2.*E0T*gb**3*lamb/(pi*m0c2)*factor_phis
        try:
            wmx  = sqrt(wmx) # T.Wangler (6.28) wmx on sepratrix
        except ValueError as ex:
            print('No energy acceptance at 1st gap! Maybe wrong frequency?')
            raise ex
        # [MeV] DW on separatrix (DE = DT a.k.a. DW)
        DWmx = wmx*m0c2       
        # Dp/p on separatrix, Dp2pmx = gamma/(gamma*gamma-1)*Dwmx
        Dp2pmx = conv.wToDp2p(wmx) # Dp/p on separatrix
        phi_1 = -phis              # [rad]
        phi_2 = 2.*phis            # [rad] Naehrung T.Wangler pp.178
        psi   = 3.*fabs(phis)      # [rad]

        # SMALL amplitude oscillations (T.Wangler pp.184)
        omgl0zuomg = sqrt(E0T*lamb*sin(-phis)/(2*pi*m0c2*gamma**3*beta))
        omgl0      = omgl0zuomg*2.*pi*freq   # [Hz]
        # Calulation of longitudinal parameters:
        #     *) Vorgabe ist delta-T/T (DT2T) kinetic energy spread
        #     *) daraus w0 als energy spread normiert auf m0c2.
        #     *) w0root ist der Wurzelaudruck in Wangler's Formel
        #     *) Dphi0 folgt aus DT2T und w0root
        #     *) emittanz folgt als Produkt von Dphi0*w0
        #     *) {Dphi,w}-space: emitw = Dphi0*w0 [rad]
        w0       = (gamma-1.)*DT2T
        w0root   = sqrt(E0T*gb**3*lamb*sin(-phis)/(2.*pi*m0c2))
        Dphi0    = (gamma-1.)/w0root*DT2T # delta-phase-intersect
        emitw    = Dphi0*w0
        betaw    = emitw/w0**2            # beta twiss
        gammaw   = 1./betaw               # gamma twiss

        # longitudinal acceptance check (always done)
        if wmx <= w0:
            si,sm,sf = node.position
            warnings.showwarning(
                'out of energy acceptance @ s={:.1f} [m]'.format(si),
                UserWarning,'setutil.py',
                'waccept()')

        # {z-dp/p}-space
        z0,Dp2p0,emitz,betaz = conv.wtoz((Dphi0,w0,emitw,betaw))
        gammaz = 1./betaz
    # TODO: use y2,y3 from Twiss or not needed? or is acceptance correct like this!
        Dp2pAcceptance = Dp2pmx
        zAcceptance    = abs(conv.DphiToz(psi))
        res =  dict(
                # {Dphi,w}
                emitw    = emitw,       # emittance{Dphi,w} [rad]
                betaw    = betaw,       # beta twiss [rad]
                gammaw   = gammaw,      # gamma twiss [1/rad]
                wmx      = wmx,         # separatrix: max w
                w0       = w0,          # w-Axenabschnitt
                Dphi0    = Dphi0,       # ellipse dphi-int (1/2 axis)
                phi_1    = phi_1,       # separatrix: max pos. phase
                phi_2    = phi_2,       # separatrix: max neg. phase
                psi      = psi,         # separatrix: bunch length [rad]
                omgl0    = omgl0,       # synchrotron oscillation [Hz]
                # {z,Dp2p}
                betaz    = betaz,       # twiss beta [m/rad]
                gammaz   = gammaz,      # twiss gamma [1/m]
                emitz    = emitz,       # emittance in {z,dp/p} space [m*rad]
                Dp2pmx   = Dp2pmx,      # max D/p on separatrix
                Dp2p0    = Dp2p0,       # ellipse dp/p-int (1/2 axis)
                z0       = z0,          # ellipse z-int    (1/2 axis) [m]
                Dp2pAcceptance = Dp2pAcceptance,
                zAcceptance    = zAcceptance,
                # {Dphi,DW}
                DWmx      = DWmx)       # separatrix: max W in [MeV]

        PARAMS['emitw']   = emitw
        PARAMS['emitz']   = emitz
        PARAMS['betaw']   = betaw
        PARAMS['betaz']   = betaz
        PARAMS['alfaz']   = 0.0
        PARAMS['DWmx']    = DWmx
        PARAMS['wmx']     = wmx
        PARAMS['psi']     = psi
        PARAMS['phi_2']   = phi_2
        PARAMS['phi_1']   = phi_1
        PARAMS['Dp2p0']   = Dp2p0
        PARAMS['Dphi0']   = Dphi0
        PARAMS['z0']      = z0
        PARAMS['omgl0']   = omgl0
        PARAMS['Dp2pmx']  = Dp2pmx
        PARAMS['w0']      = w0
        PARAMS['Dp2pAcceptance'] = Dp2pAcceptance
        PARAMS['zAcceptance']    = zAcceptance
        
        # now we can calculate the Twiss objects at injection
        alfaw = 0. # always for longitudinal
        twx = Twiss(PARAMS['betax_i'], PARAMS['alfax_i'], PARAMS['emitx_i'])
        twy = Twiss(PARAMS['betay_i'], PARAMS['alfay_i'], PARAMS['emity_i'])
        tww = Twiss(PARAMS['betaw'], alfaw, PARAMS['emitw'])
        twz = Twiss(PARAMS['betaz'],alfaw,PARAMS['emitz'])
        PARAMS['twiss_x_i'] = twx
        PARAMS['twiss_y_i'] = twy
        PARAMS['twiss_w_i'] = tww
        PARAMS['twiss_z_i'] = twz

    else:
        # assume no acceleration
        FLAGS['dWf'] = 0
        res = {}
        # we can calculate the Twiss objects at injection ...
        twx = Twiss(PARAMS['betax_i'], PARAMS['alfax_i'], PARAMS['emitx_i'])
        twy = Twiss(PARAMS['betay_i'], PARAMS['alfay_i'], PARAMS['emity_i'])
        # ... and use dummies
        tww = twz = Twiss(1.,0.,1.)
        PARAMS['twiss_x_i'] = twx
        PARAMS['twiss_y_i'] = twy
        PARAMS['twiss_w_i'] = tww
        PARAMS['twiss_z_i'] = twz
        
    PARAMS['twiss_x_i()'] = twx()
    PARAMS['twiss_xyi()'] = twy()
    PARAMS['twiss_w_i()'] = tww()
    PARAMS['twiss_z_i()'] = twz()

    return res
    
#TODO: integrate sigmas into Twiss
def sigmas(alfa,beta,epsi):
    """ calculates sigmas from twiss-alpha, -beta and -emittance """
    gamma  = (1.+ alfa**2)/beta
    sigma  = sqrt(epsi*beta)
    sigmap = sqrt(epsi*gamma)
    return sigma,sigmap

def show_data_from_elements():
    sectionIDs = LATTICE      # lATTICE has now the list of sectionIDs
    eIDsps = SECTIONS['uniqueIDs']
    DEBUG_OFF(sectionIDs)
    DEBUG_OFF(eIDsps)
    DEBUG_OFF(ELEMENTS)
    types = ['QF','QD','QFth','QDth','QFthx','QDthx','RFG','RFC']
    for sec in sectionIDs:
        elementIDs = eIDsps[sec]        
        for type in types:
            for elementID in elementIDs:
                element = ELEMENTS[elementID]
                if type == element['type']:
                    dictprnt(element,text=elementID,end='')

def collect_data_for_summary(lattice):
    if True:
        SUMMARY['use emittance growth']            =  FLAGS['egf']
        SUMMARY['use sigma tracking']              =  FLAGS['sigma']
        SUMMARY['use emittance growth']            =  FLAGS['egf']
        SUMMARY['use ring lattice']                =  FLAGS['periodic']
        SUMMARY['use express']                     =  FLAGS['express']
        SUMMARY['use aperture']                    =  FLAGS['useaper']
        SUMMARY['accON']                           =  False if  FLAGS['dWf'] == 0 else  True
        SUMMARY['lattice version']                 =  PARAMS['lattice_version']
        SUMMARY['(N)sigma']                        =  PARAMS['nbsigma']
        SUMMARY['injection energy [MeV]']          =  PARAMS['injection_energy']
        SUMMARY['(sigx )i*   [mm]']                =  PARAMS['twiss_x_i'].sigmaH()*1.e3
        SUMMARY["(sigx')i* [mrad]"]                =  PARAMS['twiss_x_i'].sigmaV()*1.e3
        SUMMARY['(sigy )i*   [mm]']                =  PARAMS['twiss_y_i'].sigmaH()*1.e3
        SUMMARY["(sigy')i* [mrad]"]                =  PARAMS['twiss_y_i'].sigmaV()*1.e3
        SUMMARY["emit{x-x'}[mrad*mm]"]             =  PARAMS['emitx_i']*1.e6
        SUMMARY["emit{y-y'}[mrad*mm]"]             =  PARAMS['emity_i']*1.e6
        SUMMARY['(delta-T/T)i spread']             =  '{:8.2e} kinetic energy'.format(PARAMS['DT2T'])
    
    if FLAGS['dWf'] == 1:
        SUMMARY['separatrix: DW-max*[MeV]']        =  '{:8.2e} energy'.format(PARAMS['DWmx'])
        SUMMARY['separatrix: w-max*   [%]']        =  '{:8.2e} delta-gamma'.format(PARAMS['wmx']*1.e2)
        SUMMARY['separatrix: Dphi*  [deg]']        =  '{:8.2f}, {:6.2f} to {:6.2f}'.format(degrees(PARAMS['psi']),degrees(PARAMS['phi_2']),degrees(PARAMS['phi_1']))
        SUMMARY['separatrix: Dp/p-max [%]']        =  '{:8.2e} impulse'.format(PARAMS['Dp2pmx']*100.)
        SUMMARY['emit{z-Dp/p}*  [mm]']             =  '{:8.2e}'.format(PARAMS['emitz']*1.e3)
        SUMMARY['emit{phi-w}*  [rad]']             =  '{:8.2e}'.format(PARAMS['emitw'])
        SUMMARY['(Dp/p)i spread*']                 =  '{:8.2e} impulse'.format(PARAMS['Dp2p0'])
        SUMMARY['(phi)i spread* [rad]']            =  '{:8.2e} phase'.format(PARAMS['Dphi0'])
        SUMMARY['(z)i spread*    [m]']             =  '{:8.2e} bunch'.format(abs(PARAMS['z0']))
        SUMMARY['sync.oscillation* [MHz]']         =  PARAMS['omgl0']*1.e-6
        SUMMARY['(w)i spread']                     =  '{:8.2e} delta-gamma, dE/E0'.format(PARAMS['w0'])
    else:
        SUMMARY['separatrix:']                     =  '{}'.format('NO acceleration')
    return

def I0(x):
    """
    Modified Bessel function I of integer order 0
    ref.: Hanbook of Mathematical Functions, M.Abramowitz & I.A.Stegun
    """
    t = x/3.75
    if 0. <= x and x < 3.75:
        t2 = t*t
        res = 1.
        res+= 3.5156229*t2
        res+= 3.0899424*t2*t2
        res+= 1.2067492*t2*t2*t2
        res+= 0.2659732*t2*t2*t2*t2
        res+= 0.0360768*t2*t2*t2*t2*t2
        res+= 0.0045813*t2*t2*t2*t2*t2*t2
        DEB.get('OFF')('(I0,x )=({},{})'.format(res,x))
    elif 3.75 <= x:
        tm1 = 1./t
        res = 0.39894228
        res+= 0.01328529*tm1
        res+= 0.00225319*tm1*tm1
        res-= 0.00157565*tm1*tm1*tm1
        res+= 0.00916281*tm1*tm1*tm1*tm1
        res-= 0.02057706*tm1*tm1*tm1*tm1*tm1
        res+= 0.02635537*tm1*tm1*tm1*tm1*tm1*tm1
        res-= 0.01647633*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        res+= 0.00392377*tm1*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        try:
            res = res*exp(x)/sqrt(x)
        except OverflowError as ex:
            print('Bessel-function I0 overflow: (arg = {:6.3f})'.format(x))
            # sys.exit(1)
            raise ex
            
    return res

def I1(x):
    """
    Modified Bessel function I of integer order 1
    ref.: Hanbook of Mathematical Functions, M.Abramowitz & I.A.Stegun
    """
    t = x/3.75
    if 0. <= x and x < 3.75:
        t2 = t*t
        res = 0.5
        res+= 0.87890059*t2
        res+= 0.51498869*t2*t2
        res+= 0.15084934*t2*t2*t2
        res+= 0.02658733*t2*t2*t2*t2
        res+= 0.00301532*t2*t2*t2*t2*t2
        res+= 0.00032411*t2*t2*t2*t2*t2*t2
        res = res*x
        DEB.get('OFF')('(I1,x )=({},{})'.format(res,x))
    elif 3.75 <= x:
        tm1 = 1/t
        res = 0.39894228
        res-= 0.03988024*tm1
        res-= 0.00362018*tm1*tm1
        res+= 0.00163801*tm1*tm1*tm1
        res-= 0.01031555*tm1*tm1*tm1*tm1
        res+= 0.02282967*tm1*tm1*tm1*tm1*tm1
        res-= 0.02895312*tm1*tm1*tm1*tm1*tm1*tm1
        res+= 0.01787654*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        res-= 0.00420059*tm1*tm1*tm1*tm1*tm1*tm1*tm1*tm1
        try:
            res = res*exp(x)/sqrt(x)
            DEB.get('OFF')('(I1,x )=({},{})'.format(res,x))
        except OverflowError as ex:
            print('Bessel-function I1 overflow: (arg = {6.3f})'.format(x))
            # sys.exit(1)
            raise ex
    return res

def ellicp(xy,alfa,beta,emit):
    """ convert twiss parameters to plot parameters """
    gamma = (1.+alfa**2)/beta
    tilt = degrees(0.5*atan(2*alfa/(gamma-beta)))  # see CERN's Formelsammlung
    a = sqrt(emit*beta)
    b = sqrt(emit/beta)
    # return matplot.patches.Ellipse(origin=xy,width=a,height=b,angle=tilt) arguments
    return (xy,a,b,tilt)

# TODO move marker actions to own module
def elli_sxy_action(node,on_injection=False):
    """ display x- and y-phase-space ellipses """
    if on_injection:
        s = 0.0
        ax = PARAMS['alfax_i']
        bx = PARAMS['betax_i']
        ay = PARAMS['alfay_i']
        by = PARAMS['betay_i']

    else:
        twiss,s = node['twiss']

        ax = twiss[Ktw.ax]
        bx = twiss[Ktw.bx]
        ay = twiss[Ktw.ay]
        by = twiss[Ktw.by]

    org = (0,0)
    ellix = ellicp(org,ax,bx,PARAMS['emitx_i'])
    elliy = ellicp(org,ay,by,PARAMS['emity_i'])

    ells = [Ellipse(*ellix,color='blue',fill=False),Ellipse(*elliy,color='red',fill=False)]

    fig, ax = plt.subplots()
    fig.suptitle('phase-space {{[m],[rad]}} @ s={:6.2f}[m]'.format(s))
    fig.legend(ells,("{x,x'}","{y,y'}"),loc=1)

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)

    x1 = sqrt(PARAMS['emitx_i']*PARAMS['betax_i'])
    x2 = sqrt(PARAMS['emity_i']*PARAMS['betay_i'])
    xmax = max(x1,x2)
    gammax = (1.+PARAMS['alfax_i']**2)/PARAMS['betax_i']
    gammay = (1.+PARAMS['alfay_i']**2)/PARAMS['betay_i']
    y1 = sqrt(PARAMS['emitx_i']*gammax)
    y2 = sqrt(PARAMS['emity_i']*gammay)
    ymax = max(y1,y2)
    # scale = 0.6
    scale = 2.0
    plt.xlim(-xmax*scale, xmax*scale)
    plt.ylim(-ymax*scale, ymax*scale)

#TODO: more marker-actions
# def sigma_x_action(*args):
#     # DEBUG_MODULE('(sigma)x @ z {:8.4f}[m] = {:8.4f}[mm]'.format(KEEP['z'],KEEP['sigma_x']*1.e3))
#     SUMMARY['z {:8.4f}[m] sigma-x [mm]'.format(KEEP['z'])] = KEEP['sigma_x']*1.e3
#     PARAMS['sigma-x({:0=6.2f})'.format(KEEP['z'])] = KEEP['sigma_x']*1.e3
# def sigma_y_action(*args):
#     # DEBUG_MODULE('(sigma)y @ z {:8.4f}[m] = {:8.4f}[mm]'.format(KEEP['z'],KEEP['sigma_y']*1.e3))
#     SUMMARY['z {:8.4f}[m] sigma-y [mm]'.format(KEEP['z'])] = KEEP['sigma_y']*1.e3
#     PARAMS['sigma-y({:0=6.2f})'.format(KEEP['z'])] = KEEP['sigma_y']*1.e3
# def tkin_action(*args):
#     SUMMARY['z {:8.4f}[m]   Tkin [MeV]'.format(KEEP['z'])] = KEEP['Tkin']
#     PARAMS['Tkin({:0=6.2f})'.format(KEEP['z'])] = KEEP['Tkin']
# 
# """
#  (global) MRKR_ACTIONS: dictionary of possible actions attached to a marker
# """
# MRKR_ACTIONS = dict(
#             sigma_x     = sigma_x_action,
#             sigma_y     = sigma_y_action,
#             Tkin        = tkin_action,
#             show_elli   = elli_sxy_action,
#             )

# utilities
def k0prot(gradient=0.,tkin=0.):
    """
    Quadrupole strength as function of kin. energy and gradient (only for protons!)
    IN:
        gradient: in [T/m]
        tkin: kinetic energy in [MeV]
    OUT:
        k in [1/m^2]
    """
    if (tkin >= 0.):
        prot = Proton(tkin)
        k    = gradient/prot.brho
        return k
    else:
        print('k0prot() called with negative kinetic energy? - STOP')
        sys.exit(1)

def scalek0prot(k0=0.,tki=0.,tkf=0.):
    """Scale Quadrupole strength k0 for increase of kin. energy from tki to tkf  (only for protons!)"""
    bgi  = Proton(tki).gamma_beta
    bgf  = Proton(tkf).gamma_beta
    k= k0 * bgi/bgf
    return k

def dBdxprot(k0=0.,tkin=0.):
    """
    B-field gradient from quadrupole gradient for given quadrupole strength k0 and kin. energy tkin (only for protons!)
    IN:
        k0 qudrupole strength [1/m^2]
        tkin kin.energy [MeV]
    OUT:
        dB/dz in [T/m]
    """
    if (tkin >= 0.):
        return k0 * Proton(tkin).brho
    else:
        print('dBdxprot() called with negative kinetic energy? - STOP')
        sys.exit(1)

def objprnt (what,text='',filter=[]):
    """Custom helper to print objects as dictionary"""
    template = '============================================'
    lt  = len(template)
    lx  = len(text)
    p1 = int((lt-lx)/2)
    p2 = int((lt+lx)/2)
    if p1 < 0:
        ueberschrift = text
    else:
        ueberschrift = template[:p1]+text+template[p2:]
    print('\n          '+ueberschrift)
    for k,v in sorted(what.__dict__.items()):
        if k in filter:
            continue
        if k == 'matrix':
            print(k.rjust(30),':\n',v)
        else:
            print(k.rjust(30),':',v)
    return

def dictstring(what,text='',filter=[],njust=35):
    """Custom helper to print dictionaries (clever!?)"""
    def asString(v):
        txt = ''
        fmt0  = '{:8.6g} '
        fmt1 = '{} '
        fmt2 = '{:s} '
        if isinstance(v,bool):
            txt += fmt1.format(v)
        elif isinstance(v,str):
            txt += fmt2.format(v)
        elif isinstance(v,float) or isinstance(v,int) or isinstance(v,complex):
            txt += fmt0.format(v)
        else:
            txt += fmt1.format(v)
        return txt

    template = '=============================================='
    res = ''
    lt  = len(template)
    lx  = len(text)
    p1 = int((lt-lx)/2)
    p2 = int((lt+lx)/2)
    if p1 < 0:
        ueberschrift = text
    else:
        ueberschrift = template[:p1]+' {} '.format(text)+template[p2:]
    # print('            '+ueberschrift)
    res+= '            {}\n'.format(ueberschrift)

    fmt  = '{:>'+'{}s'.format(njust)+'} : '
    for k,v in sorted(what.items()):
        if k in filter:
            continue
        vars = ''
        if isinstance(v,tuple):
            for i in v: vars += asString(i)
        else:
            vars += asString(v)
        # print(fmt.format(k)+vars)
        res+=fmt.format(k)+'{}\n'.format(vars)
    return res

def dictprnt(what,text='',filter=[],njust=35,end='\n'):
    print(dictstring(what,text,filter,njust),end=end)

def printv(level,*args):
    """Multilevel printing with verbose flag"""
    verbose = FLAGS['verbose']
    if verbose >= level and not FLAGS['KVout']:
        print(*args)

def tblprnt(headr,records):
    """
    Custom helper to print a nice table to memory (clever!?)
    IN:
        headr   = table header [...]
        records = table rows [[...]]
    OUT:
        s = the table as a string
    """
    rows = []; s=''
    rows.append(headr)
    for record in records:
        row = record
        rows.append(row)
    widths = [max(map(len, col)) for col in zip(*rows)]
    for row in rows:
            s+=" | ".join((val.ljust(width) for val,width in zip(row, widths)))+'\n'
    return s

def mxprnt(matrix):
    """Simple matrix print"""
    s = [['{:+.3e}  '.format(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = ''.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    return '\n'.join(table)

def arrprnt(array,fmt='{:8.4f}, ',txt=''):
    """Simple array print"""
    print(txt,end='')
    for val in array:
        print(fmt.format(val),end='')
    print('')

def wille():
    return {
        'k_quad_f':1.2,
        'length_quad_f':0.2,
        'k_quad_d':1.2,
        'length_quad_d':0.4,
        'bending_radius':3.8197,
        # 'bending_radius':1.e4,
        'dipole_length':1.5,
        'drift_length':0.55
    }

def test0():
    print('--------------------------Test0---')
    dictprnt( PARAMS,text=' PARAMS')

    print('Sollteilchen\n'+ PARAMS['sollteilchen'].string())
    print('Proton(tkin=5.)\n'+Proton(tkin=5.).string())
    print('Electron(tkin=5.)\n'+Electron(tkin=5.).string())

    dictprnt(wille(),text='wille')
    kqf = wille()['k_quad_f']
    tk  = 5.
    headr = ['kq[1/m^2]','tk[Mev]','dBdxprot(kqf,tk)[T/m]']
    records = [['{:4.4f}'.format(kqf),'{:4.4f}'.format(tk),'{:4.4f}'.format(dBdxprot(k0=kqf,tkin=tk))]]
    print('\n'+tblprnt(headr,records))

def test1():
    print('--------------------------Test1---')
    args = ellicp((0,0),0.5,100.,1.e-6)
    ells = [Ellipse(*args,fill=False,color='red')]

    # fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    fig, ax = plt.subplots()
    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)

    width  = args[1]
    height = args[2]
    scale = 0.7
    ax.set_xlim(-width*scale, width*scale)
    ax.set_ylim(-height*scale, height*scale)
    plt.show()

def test2():
    print('--------------------------Test2---')
    print('test particle energy adjustment...')
    p = PARAMS['sollteilchen']
    print(repr(p)+':')
    print(p.string())
    p1 = p(100.)
    print(repr(p1)+':')
    print(p1.string())
    try:
        p2 = p(-10.)
    except ValueError:
        print("This was a test: Particle energy deliberately set negative!")

def test3():
    print('--------------------------Test3---')
    names = ('t','SIN','COS')
    fn = Functions(names)
    for ph in range(0,361):
        phr = radians(ph)
        fn.append(ph,(sin(phr),cos(3*phr)))
    
    t= [fn(i,'t')   for i in range(fn.nbpoints)]
    S= [fn(i,'SIN') for i in range(fn.nbpoints)]
    C= [fn(i,'COS') for i in range(fn.nbpoints)]
    
    for i in range(fn.nbpoints):
        print(t[i],S[i],C[i])
    
    plt.plot(S,C)
    plt.show()

if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()

