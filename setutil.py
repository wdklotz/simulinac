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

import sys,traceback
import scipy.constants as C
from math import pi,sqrt,sin,cos,radians,degrees,fabs,exp,atan
import logging, pprint
from enum import IntEnum
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import warnings

# DEBUG
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

def DEBUG(string,arg='',end='\n'):
    """
    Print debug message
    IN:
        string: text to print
        arg:  values to print
    """
    if isinstance(arg,list):
        # print('DEBUG: {} \nlist={}'.format(string,arg))
        pp   = pprint.PrettyPrinter(indent=4)  # use pprint module
        sarg = pp.pformat(arg)

        print('DEBUG: {} typ(list) {}'.format(string,sarg),end=end)
    elif isinstance(arg,dict):
        # print('DEBUG: {} \ndict={}'.format(string,arg))
        pp   = pprint.PrettyPrinter(indent=4,width=60)  # use pprint module
        sarg = pp.pformat(arg)
        print('DEBUG: {} typ(dict) {}'.format(string,sarg),end=end)
    else:
        print('DEBUG: {}{}'.format(string,arg),end=end)

# Logger
ch        = logging.StreamHandler()     # console handler
ch.setLevel(logging.DEBUG)              # set handler level
formatter = \
    logging.Formatter("%(levelname)s: %(filename)s[%(lineno)d] %(message)s")
ch.setFormatter(formatter)              # set handler's format
logger    = logging.getLogger("logger")
logger.addHandler(ch)                   # add handler to logger

# x        x'        y        y'        z       z'=dp/p   T        dT        S        dS
XKOO = 0;XPKOO = 1;YKOO = 2;YPKOO = 3;ZKOO = 4;ZPKOO = 5;EKOO = 6;DEKOO = 7;SKOO = 8;LKOO = 9
class K(IntEnum):                                # enum.IntEnum since Python 3.4
    """ Koordinaten for track points (1x10)"""
    x  = XKOO
    xp = XPKOO
    y  = YKOO
    yp = YPKOO
    z  = ZKOO
    zp = ZPKOO
    T  = EKOO
    dT = DEKOO
    S  = SKOO
    dS = LKOO

class K6(IntEnum):
    """ Koordinaten for twiss funtions (1x4) """
    bx = 0      # twiss-beta
    ax = 1      # twiss-alpha
    gx = 2      # twiss-gamma
    by = 3
    ay = 4
    gy = 5
    s  = 6      # abzisse for twiss functions

# DEFAULTS "FLAGS" & "PARAMS"
FLAGS  = dict(
        periodic             = False,            # periodic lattice? default
        egf                  = False,            # emittance grow flag default
        sigma                = True,             # beam sizes by sigma-tracking
        KVout                = False,            # print a dictionary of Key-Value pairs, no display
        dWf                  = 1.,               # acceleration on/off flag 1=on,0=off
        verbose              = 0,                # print flag default = 0
        express              = False,            # use express version of thin quads
        useaper              = False,            # use aperture check for quads and rf-gaps
        bucket               = False,            # plot bucket
        csTrak               = True,             # plot CS trajectories
        pspace               = False             # plot CS twiss ellipses at entrance
        )
PARAMS = dict(
        lichtgeschwindigkeit = 299792458.,       # [m/s] const
        elementarladung      = 1.602176565e-19,  # [coulomb] const
        proton_mass          = 938.272,          # [MeV/c**2] const
        electron_mass        = 0.5109989,        # [MeV/c**2] const
        EzAvg                = 2.87,             # [MV/m] default average E-field on axis
        gap                  = 0.022,            # [m] default
        cavity_laenge        = 0.08,             # [m] default
        phisoll              = -30.,             # [deg] default
        frequenz             = 816.e6,           # [Hz] default
        injection_energy     = 50.,              # [MeV] default
        qf_gradient          = 16.0,             # [T/m] default
        qd_gradient          = 16.0,             # [T/m] default
        quad_bore            = 0.02,             # [m] Vorgabe quadrupole bore radius
        aperture             = None,             # default aperture = no aperture
        n_coil               = 30,               # nbof coil windings
        n_sigma              = 3,                # nboff beam sigmas to stay clear of aperture
        emitx_i              = 2.0e-6,           # [m*rad] Vorgabe emittance entrance
        emity_i              = 2.0e-6,           # [m*rad] Vorgabe emittance entrance
        emitw_i              = 0.2e-6,           # [rad] Vorgabe emittance entrance
        betax_i              = 2.800,            # [m] Vorgabe twiss beta entrance
        betay_i              = 0.200,            # [m] Vorgabe twiss beta entrance
        alfax_i              = 0.0,              # Vorgabe twiss alpha entrance
        alfay_i              = 0.0,              # Vorgabe twiss alpha entrance
        alfaw_i              = 0.0,              # Vorgabe twiss alpha entrance
        nbof_slices          = 10,               # default number of slices
        mapset               = frozenset(['t3d','simple','base','ttf','dyn']), #gap-models
        )
"""
 (global) KEEP: dict to keep tracking results
"""
KEEP = dict(z=0.,sigma_x=0.,sigma_y=0.,Tkin=0.)

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
        self.v          = self.beta * PARAMS['lichtgeschwindigkeit']    # velocity [m/s]
        self.brho       = 1.e+6 / PARAMS['lichtgeschwindigkeit'] * self.gamma_beta * self.e0 # [T*m]
        self.name       = name
        self.m0c2       = self.e0
        self.m0c3       = self.e0 * PARAMS['lichtgeschwindigkeit']
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
    def trtf(self,gap,fRF):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*fRF*gap / (self.beta* PARAMS['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def __call__(self,tkin):  # make callable to change energy
        self._set_self(tkin=tkin,mass=self.e0,name=self.name)
        return self

class Proton(Particle):
    def __init__(self,tkin= PARAMS['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass= PARAMS['proton_mass'],name='proton')

class Electron(Particle):
    def __init__(self,tkin= PARAMS['injection_energy']):
        super(Electron,self).__init__(tkin=tkin,mass= PARAMS['electron_mass'],name='electron')

# Sollteichen
PARAMS['sollteilchen'] = Proton()
PARAMS['wellenlänge'] = PARAMS['lichtgeschwindigkeit']/PARAMS['frequenz']

class WConverter(object):
    """
    Converter to switch between different longitudinal phase space coordinates
    """
    def __init__(self,tk,freq=PARAMS['frequenz']):
        self.pi             = C.pi
        self.lamb           = C.c/freq           # [m]
        self.m0c2,unit,prec = C.physical_constants['proton mass energy equivalent in MeV']
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

    def emitwToemitz(self,emitw):
        """ emittance[w(x)Dphi] [rad] to emittance[z(x)Dp2p] [m] """
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

## Long. Emittance
def waccept(node):
    """
    Central to calculate longitudinal phase space ellipse parameters nach T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0*Dphi0 = ellipse_area/pi
    IN
        node: the 1st rf-gap at the linac entrance
    """
    if node is not None and FLAGS['dWf']:
        emitw_i   = PARAMS['emitw_i']    # [rad]
        E0T       = node.EzAvg*node.tr   # [MV/m]
        particle  = node.particle
        phis      = node.phis            # [rad]
        lamb      = node.lamb            # [m]
        freq      = node.freq            # [Hz]
        m0c2      = particle.e0          # [MeV]
        gb        = particle.gamma_beta
        beta      = particle.beta
        gamma     = particle.gamma
        tkin      = particle.tkin
        conv      = WConverter(tkin)

        # large amplitude oscillations (T.Wangler pp. 175)
        # NOTE: w == Dgamma == normalized energy deviation
        factor_phis = phis*cos(phis)-sin(phis)
        wmx  = sqrt(2.*E0T*gb**3*lamb/(pi*m0c2)*factor_phis) # T.Wangler (6.28) wmx on sepratrix
        DWmx = wmx*m0c2       # [MeV] DW on separatrix (DE = DT == DW)
        # NOTE: Dp2pmx = gamma/(gamma*gamma-1)*wmx # Dp/p on separatrix
        Dp2pmx = conv.wToDp2p(wmx) # Dp/p on separatrix
        phi_1 = -phis              # [rad]
        phi_2 = 2.*phis            # [rad] Naehrung T.Wangler pp.178
        psi   = 3.*fabs(phis)      # [rad]

        # small amplitude oscillations (T.Wangler pp.184)
        omgl0zuomg = sqrt(E0T*lamb*sin(-phis)/(2*pi*m0c2*gamma**3*beta))
        omgl0      = omgl0zuomg*2.*pi*freq   # [Hz]
        # {Dphi,w}-space   
        # NOTE: emitw_i = Dphi0*w0 = w0root*Dphi0**2  [rad]
        w0root   = sqrt(E0T*gb**3*lamb*sin(-phis)/(2.*pi*m0c2))
        Dphi0    = sqrt(emitw_i/w0root)     # delta-phase-intersect
        w0       = w0root*Dphi0             # w-intersect 
        betaw    = emitw_i/w0**2            # beta twiss
        gammaw   = 1./betaw                 # gamma twiss

        # longitudinal acceptance check (always done)
        if wmx <= w0:
            si,sm,sf = node.position
            warnings.showwarning(
                'out of energy acceptance @ s={:.1f} [m]'.format(si),
                UserWarning,'setutil.py',
                'waccept()')

        # {z,dp/p}-space
        # z0       = -Dphi0*beta*lamb/(2.*pi)        # z [m]
        # Dp2p0    = w0/(beta*beta*gamma)            # delta-p/p
        # emitz    = -z0*Dp2p0                       # emittance {z,Dp/p} [m]
        # betaz    = emitz/Dp2p0**2                  # beta [m]
        z0,Dp2p0,emitz,betaz = conv.wtoz((Dphi0,w0,emitw_i,betaw))
        gammaz = 1./betaz
#todo: use y2,y3 from Twiss or not needed? Acceptance is correct like this!
        Dp2pAcceptance = Dp2pmx
        zAcceptance    = abs(conv.DphiToz(psi))
        res =  dict(
                # {Dphi(x)w}
                betaw    = betaw,       # beta twiss [rad]
                gammaw   = gammaw,      # gamma twiss [1/rad]
                wmx      = wmx,         # separatrix: max w
                w0       = w0,          # w-Axenabschnitt
                Dphi0    = Dphi0,       # ellipse dphi-int (1/2 axis)
                phi_1    = phi_1,       # separatrix: max pos. phase
                phi_2    = phi_2,       # separatrix: max neg. phase
                psi      = psi,         # separatrix: bunch length [rad]
                omgl0    = omgl0,       # synchrotron oscillation [Hz]
                # {z(x)Dp2p}
                betaz    = betaz,       # twiss beta [m]
                gammaz   = gammaz,      # twiss gamma [1/m]
                emitz    = emitz,       # emittance in {z,dp/p} space [m]
                Dp2pmx   = Dp2pmx,      # max D/p on separatrix
                Dp2p0    = Dp2p0,       # ellipse dp/p-int (1/2 axis)
                z0       = z0,          # ellipse z-int    (1/2 axis) [m]
                Dp2pAcceptance = Dp2pAcceptance,
                zAcceptance    = zAcceptance,
                # {Dphi(x)DW}
                DWmx      = DWmx)       # separatrix: max W in [MeV]

    PARAMS['emitz']   = emitz
    PARAMS['betaw']   = betaw
    PARAMS['betaz']   = betaz
    PARAMS['DWmx']    = DWmx
    PARAMS['wmx']     = wmx
    PARAMS['psi']     = psi
    PARAMS['phi_2']   = phi_2
    PARAMS['phi_1']   = phi_1
    PARAMS['Dp2p0']   = Dp2p0
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
    tww = Twiss(PARAMS['betaw'], alfaw, PARAMS['emitw_i'])
    twz = Twiss(PARAMS['betaz'],alfaw,PARAMS['emitz'])
    PARAMS['twiss_x_i'] = twx
    PARAMS['twiss_y_i'] = twy
    PARAMS['twiss_w_i'] = tww
    PARAMS['twiss_z_i'] = twz

    return res
#todo: integrate sigmas into Twiss
def sigmas(alfa,beta,epsi):
    """ calculates sigmas from twiss-alpha, -beta and -emittance """
    gamma  = (1.+ alfa**2)/beta
    sigma  = sqrt(epsi*beta)
    sigmap = sqrt(epsi*gamma)
    return sigma,sigmap

# (global) SUMMARY: dictionary for summary
SUMMARY = {}

def collect_data_for_summary(lattice):
    class Filter:
        def __init__(self,func):
            self.func = func
        def __ror__(self,iterable):     # INFIX "|" operator
            for obj in iterable:
                if self.func(obj):
                    yield obj
    class Apply:
        def __init__(self,func):
            self.func = func
        def __ror__(self,iterable):     # INFIX "|" operator
            return self.func(iterable)

    def elements_in_lattice():
        '''
        Filter elements of class <typ> and section <sec> from lattice
        IN:
            lattice = object        [Lattice]
            typ     = element class [string]
            sec     = section name  [string]
        OUT:
            list of filtered elements

        NOTE: this functional implementation is taken from:
            https://code.activestate.com/recipes/580625-collection-pipeline-in-python/
        '''
        def predicate(element):
            try:
                test = (type(element).__name__ == typ and element.section == sec)
            except AttributeError:
                test = (type(element).__name__ == typ)  # no section tag? take all!
            return test

        # NOTE: here I use the INFIX operator '|' like a UNIX pipe
        List     = Apply(list)
        Selector = Filter(predicate)
        return lattice.seq | Selector | List

    def elements_in_section():
        """Remove duplicate elements of same type in a section"""
        elements = elements_in_lattice()
        new_elements = []
        seen = set()             # helper to eliminate duplicate entries
        for itm in elements:
            label = itm.label
            if label in seen:
                continue
            else:
                seen.add(label)
                new_elements.append(itm)
        return new_elements

    # body
    sections =  PARAMS['sections']                   # comes from INPUT
    if not FLAGS['sections']: sections = ['*']       # section wildcart
    types = ['QF','QD','QFth','QDth','QFthx','QDthx']
    for sec in sections:
        for typ in types:
            elements = elements_in_section()
            for itm in elements:
                k0       = itm.k0
                dBdz     = k0*itm.particle.brho
                length   = itm.length
                aperture = itm.aperture
                # SUMMARY['{2} [{1}.{0}]    k0 [m^-2]'.format(sec,typ,itm.label)] = k0
                SUMMARY['{2} [{1}.{0}]    dBdz[T/m]'.format(sec,typ,itm.label)] = dBdz
                SUMMARY['{2} [{1}.{0}]       B0*[T]'.format(sec,typ,itm.label)] = dBdz*aperture
                SUMMARY['{2} [{1}.{0}]    length[m]'.format(sec,typ,itm.label)] = length
                SUMMARY['{2} [{1}.{0}]  aperture[m]'.format(sec,typ,itm.label)] = aperture

    types = ['RFG']
    for sec in sections:
        for typ in types:
            elements = elements_in_section()
            for itm in elements:
                gap      = itm.gap
                EzAvg    = itm.EzAvg
                PhiSoll  = degrees(itm.phis)
                mapping  = itm.mapping
                EzPeak   = itm.EzPeak
                aperture = itm.aperture
                SUMMARY['{2} [{1}.{0}]       gap[m]'.format(sec,typ,itm.label)] = gap
                SUMMARY['{2} [{1}.{0}]  aperture[m]'.format(sec,typ,itm.label)] = aperture
                SUMMARY['{2} [{1}.{0}]  EzAvg[MV/m]'.format(sec,typ,itm.label)] = EzAvg
                SUMMARY['{2} [{1}.{0}]    phis[deg]'.format(sec,typ,itm.label)] = PhiSoll
                SUMMARY['{2} [{1}.{0}]      mapping'.format(sec,typ,itm.label)] = mapping
                SUMMARY['{2} [{1}.{0}] EzPeak[MV/m]'.format(sec,typ,itm.label)] = EzPeak

    types = ['RFC']
    for sec in sections:
        for typ in types:
            elements = elements_in_section()
            for itm in elements:
                gap      = itm.gap
                EzAvg    = itm.EzAvg
                PhiSoll  = degrees(itm.phis)
                length   = itm.length
                aperture = itm.aperture
                SUMMARY['{2} [{1}.{0}]       gap[m]'.format(sec,typ,itm.label)] = gap
                SUMMARY['{2} [{1}.{0}]  aperture[m]'.format(sec,typ,itm.label)] = aperture
                SUMMARY['{2} [{1}.{0}]  EzAvg[MV/m]'.format(sec,typ,itm.label)] = EzAvg
                SUMMARY['{2} [{1}.{0}]    phis[deg]'.format(sec,typ,itm.label)] = PhiSoll
                SUMMARY['{2} [{1}.{0}]    length[m]'.format(sec,typ,itm.label)] = length

    SUMMARY['use sigma tracking']              =  FLAGS['sigma']
    SUMMARY['use emittance growth']            =  FLAGS['egf']
    SUMMARY['use ring lattice']                =  FLAGS['periodic']
    SUMMARY['use express']                     =  FLAGS['express']
    SUMMARY['use aperture']                    =  FLAGS['useaper']
    SUMMARY['accON']                           =  False if  FLAGS['dWf'] == 0. else  True
    SUMMARY['wavelength [cm]']                 =  PARAMS['wellenlänge']*1.e2
    SUMMARY['lattice version']                 =  PARAMS['lattice_version']
    SUMMARY['frequency [MHz]']                 =  PARAMS['frequenz']*1.e-6
    SUMMARY['(N)sigma']                        =  PARAMS['n_sigma']
    SUMMARY['injection energy [MeV]']          =  PARAMS['injection_energy']
    SUMMARY['(sigx )i*   [mm]']                =  PARAMS['twiss_x_i'].sigmaH()*1.e3
    SUMMARY["(sigx')i* [mrad]"]                =  PARAMS['twiss_x_i'].sigmaV()*1.e3
    SUMMARY['(sigy )i*   [mm]']                =  PARAMS['twiss_y_i'].sigmaH()*1.e3
    SUMMARY["(sigy')i* [mrad]"]                =  PARAMS['twiss_y_i'].sigmaV()*1.e3
    SUMMARY['separatrix: DW*   [MeV]']         =  '{:8.2e}'.format(PARAMS['DWmx'])
    SUMMARY['separatrix: w*      [%]']         =  '{:8.2e}'.format(PARAMS['wmx']*1.e2)
    SUMMARY['separatrix: Dphi* [deg]']         =  '{:8.2f}, {:6.2f} to {:6.2f}'.format(degrees(PARAMS['psi']),degrees(PARAMS['phi_2']),degrees(PARAMS['phi_1']))
    SUMMARY["emit{x,x'}[mrad*mm]"]             =  PARAMS['emitx_i']*1.e6
    SUMMARY["emit{y,y'}[mrad*mm]"]             =  PARAMS['emity_i']*1.e6
    SUMMARY['emit{Dphi,w} [mrad]']             =  '{:8.2e}'.format(PARAMS['emitw_i']*1.e3)
    SUMMARY['emit{z,Dp/p}*  [mm]']             =  '{:8.2e}'.format(PARAMS['emitz']*1.e3)
    SUMMARY['(sig-Dp/p)i* [%]']                =  '{:8.2e}'.format(PARAMS['Dp2p0']*1.e2)
    SUMMARY['(sigz)i*    [cm]']                =  '{:8.2e}'.format(PARAMS['z0']*1.e2)
    SUMMARY['sync.oscillation* [MHz]']         =  PARAMS['omgl0']*1.e-6
    SUMMARY['Dp/p-max on separatrix [%]']      =  PARAMS['Dp2pmx']*100.
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
        # DEBUG_MODULE('(I0,x )',(res,x))
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
            # DEBUG_MODULE('(I0,x )',(res,x))
        except OverflowError as ex:
            print('Bessel-function I0 overflow: (arg = {:6.3f})! - STOP'.format(x))
            sys.exit(1)
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
        # DEBUG_MODULE('(I1,x )',(res,x))
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
            # DEBUG_MODULE('(I1,x )',(res,x))
        except OverflowError as ex:
            print('Bessel-function I1 overflow: (arg = {6.3f})! - STOP'.format(x))
            sys.exit(1)
    return res

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

def ellicp(xy,alfa,beta,emit):
    """ convert twiss parameters to plot parameters """
    gamma = (1.+alfa**2)/beta
    tilt = degrees(0.5*atan(2*alfa/(gamma-beta)))  # see CERN's Formelsammlung
    a = sqrt(emit*beta)
    b = sqrt(emit/beta)
    # return matplot.patches.Ellipse(origin=xy,width=a,height=b,angle=tilt) arguments
    return (xy,a,b,tilt)

# marker actions
def elli_sxy_action(on_injection=False):
    """ display x- and y-phase-space ellipses """
    if on_injection:
        s = 0.0
        ax = PARAMS['alfax_i']
        bx = PARAMS['betax_i']
        ay = PARAMS['alfay_i']
        by = PARAMS['betay_i']

    else:
        node  = args[0]
        twiss,s = node['twiss']

        ax = twiss[K6.ax]
        bx = twiss[K6.bx]
        ay = twiss[K6.ay]
        by = twiss[K6.by]

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

def sigma_x_action(*args):
    # DEBUG_MODULE('(sigma)x @ z {:8.4f}[m] = {:8.4f}[mm]'.format(KEEP['z'],KEEP['sigma_x']*1.e3))
    SUMMARY['z {:8.4f}[m] sigma-x [mm]'.format(KEEP['z'])] = KEEP['sigma_x']*1.e3
    PARAMS['sigma-x({:0=6.2f})'.format(KEEP['z'])] = KEEP['sigma_x']*1.e3
def sigma_y_action(*args):
    # DEBUG_MODULE('(sigma)y @ z {:8.4f}[m] = {:8.4f}[mm]'.format(KEEP['z'],KEEP['sigma_y']*1.e3))
    SUMMARY['z {:8.4f}[m] sigma-y [mm]'.format(KEEP['z'])] = KEEP['sigma_y']*1.e3
    PARAMS['sigma-y({:0=6.2f})'.format(KEEP['z'])] = KEEP['sigma_y']*1.e3
def tkin_action(*args):
    SUMMARY['z {:8.4f}[m]   Tkin [MeV]'.format(KEEP['z'])] = KEEP['Tkin']
    PARAMS['Tkin({:0=6.2f})'.format(KEEP['z'])] = KEEP['Tkin']
def poincare_action(*args):
    pass

"""
 (global) ACTIONS: dictionary of possible actions attached to a marker
"""
ACTIONS = dict(
            sigma_x     = sigma_x_action,
            sigma_y     = sigma_y_action,
            Tkin        = tkin_action,
            show_elli   = elli_sxy_action,
            poincare    = poincare_action
            )

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
        print('dBdxprot() called with negative kinetic energy?' - STOP)
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

def dictprnt(what,text='',filter=[],njust=35):
    """Custom helper to print dictionaries (clever!?)"""
    def asstrng(v):
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

    template = '==============================================='
    lt  = len(template)
    lx  = len(text)
    p1 = int((lt-lx)/2)
    p2 = int((lt+lx)/2)
    if p1 < 0:
        ueberschrift = text
    else:
        ueberschrift = template[:p1]+text+template[p2:]
    print('\n          '+ueberschrift)

    fmt  = '{:>'+'{}s'.format(njust)+'} : '
    for k,v in sorted(what.items()):
        if k in filter:
            continue
        vars = ''
        if isinstance(v,tuple):
            for i in v: vars += asstrng(i)
        else:
            vars += asstrng(v)
        print(fmt.format(k)+vars)
    return

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
    p2 = p(-10.)

if __name__ == '__main__':
    test0()
    test1()
    test2()

