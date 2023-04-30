#!/Users/klotz/anaconda3/bin/python3.6
#-*- coding: utf-8 -*-
# Python 2 and 3 print compatability
from __future__ import print_function   #TODO still used?
__version__='v11.0.0'
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
import scipy.constants as C
import scipy.special as SP
import matplotlib.pyplot as plt
import time
import lattice_parser_2 as parser
import unittest
from math import pi,sqrt,sin,cos,radians,degrees,fabs,exp,atan
from enum import IntEnum
from matplotlib.patches import Ellipse

# DEBUG utility used by all modules
import pprint, inspect
def PRINT_PRETTY(*args):
    frameinfo = inspect.stack()[1]    # caller's stack frame
    file=frameinfo[1]
    line=frameinfo[2]
    print('DEBUG {}, {} =======> '.format(file,line),end="")
    if len(args) == 0: print()
    else:
        for obj in args:
            pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
    return True
def PASS(*args):
    return False
DEB       = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

# MDIM: dimension of matrices
MDIM = 10
class Ktp(IntEnum):
    """ TrackPoint K-oordinaten Ktp """
    x  = 0     # x
    xp = 1     # x'
    y  = 2     # y
    yp = 3     # y'
    z  = 4     # z
    zp = 5     # z' = Dp/p
    T  = 6     # current kin energy T(s)
    dT = 7     # const 1
    S  = 8     # current position sf
    dS = 9     # const 1
    # nd = 10    # current node   TODO needed?
class Ktw(IntEnum):
    """ -K-oordinaten -tw-iss vector (1x10) """
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
# for compatability with elder code!
XKOO     = Ktp.x; XPKOO=Ktp.xp; YKOO=Ktp.y; YPKOO=Ktp.yp; ZKOO=Ktp.z
ZPKOO    = Ktp.zp; EKOO=Ktp.T; DEKOO=Ktp.dT; SKOO=Ktp.S; LKOO=Ktp.dS
RUN_MODE = {0:'ring with cavities', 1:'linac', 2:'ring w/o cavities', 3:'line w/o cavities'}
FLAGS    = dict(
    # lattice_generator may override some of these FLAGs
    dWf = 1
    )
PARAMS   = dict(
    # global constants
    clight               = C.c, # [m/s]
    elementarladung      = C.e, # [coulomb]
    proton_mass          = C.value('proton mass energy equivalent in MeV'),
    electron_mass        = C.value('electron mass energy equivalent in MeV'),
    map_set              = ['t3d','simple','base','ttf','dyn','oxal'],
    warnmx               = 10, # max warnings
    injection_energy     = 50. # default
    )              
ELEMENTS = {}
SUMMARY  = {}
class Twiss(object):
    def __init__(self, beta=1., alfa=0.,epsi=0.):
        self.beta  = beta
        self.alfa  = alfa
        self.epsi  = epsi
        self.gamma = (1.+self.alfa**2)/self.beta
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
    """ A particle class """
    def __init__(self,tkin,mass,name):
        self.lost       = False
        self._track     = None
        self._set_self(tkin,mass,name)
    def _set_self(self,tkin,mass,name):
        self.tkin       = tkin                     # kinetic energy [MeV]
        self.e0         = mass                     # rest mass [MeV]
        self.e          = self.e0+self.tkin        # total energy [MeV]
        self.gamma      = self.e/self.e0
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
    @property
    def track(self):
        return self._track
    @track.setter
    def track(self,value):
        self._track=value
    def toString(self):
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
    def __init__(self,tkin):
        super(Proton,self).__init__(tkin,PARAMS['proton_mass'],'proton')
        self.brho = 3.1297*self.gamma_beta
class Electron(Particle):
    def __init__(self,tkin):
        super(Electron,self).__init__(tkin,PARAMS['electron_mass'],'electron')
class WConverter(object):
    """
    Converter to switch between different longitudinal phase space coordinates
    W is kin energy (Wangler 1.2) (gamma-1)*m0c2
    w is Dgamma (Wangler 6.13)  DW=w*m0c2
    """
    #TODO extend to other particles not only protons: marbe!
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
        """ delta-phi [rad] to z [m] """
        z = -self.bl/self.twopi*Dphi  # trace-3d pp 4
        return z    # [m]
    def zToDphi(self,z):
        """ z [m] to delta-phi [rad] """
        Dphi = -self.twopi/self.bl*z
        return Dphi    # [rad]

    def wToDp2p(self,w):
        """ w=DW/m0c2(a.k.a. Dgamma) [] to Dp2p [] """
        Dp2p = self.g/(self.g2-1.)*w   # CERN Formelsammlung tab. 1.3
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
        # conversion phase-space coordinates {Dpi(x)w} to {z(x)Dp2p}
        """ 
        IN: args is tuple = (Dphi, w (== dT/m0c2), emittance-w, beta-w)
        OUT:     is tuple = (z, Dp2p, emittance-z, beta=z)
        """
        Dphi, w, emitw, betaw = args
        z = self.DphiToz(Dphi)
        Dp2p = self.wToDp2p(w)
        emitz = self.emitwToemitz(emitw)
        betaz = self.betawTobetaz(betaw)
        return (z,Dp2p,emitz,betaz)
class Functions(object):         #TODO better use pandas?
    """ A class to gather function-values (Ordinaten) over a common independent variable (Abszisse) """
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
    """ A (singleton) container for objects  (not used, not finished) """
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
        # self._t0 = t
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
class colors: # You may need to change color settings
    RED = '\033[31m'
    ENDC = '\033[m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
def sigmas(alfa,beta,epsi):     #TODO: integrate sigmas into Twiss- maybe?
    """ calculates sigmas from twiss-alpha, -beta and -emittance """
    gamma  = (1.+ alfa**2)/beta
    sigma  = sqrt(epsi*beta)
    sigmap = sqrt(epsi*gamma)
    return sigma,sigmap
def show_data_from_elements():  #TODO better get data fron lattice objects - maybe?
    eIDs = parser.parse().ELMIDs    # get field from parser results - ??
    for elementID in sorted(eIDs):      
        element = ELEMENTS[elementID]
        print('{} '.format(elementID),end='')
        dictprnt(element,'(MKSA units)',end='')
def collect_data_for_summary(lattice):
    # Unicode characters in python :https://pythonforundergradengineers.com/unicode-characters-in-python.html
    if True:
        SUMMARY['use emittance growth']            =  FLAGS['egf']
        SUMMARY['use sigma tracking']              =  FLAGS['sigma']
        SUMMARY['use emittance growth']            =  FLAGS['egf']
        SUMMARY['use ring lattice']                =  FLAGS['periodic']
        SUMMARY['use aperture']                    =  FLAGS['useaper']
        SUMMARY['accON']                           =  FLAGS['accON']
        SUMMARY['lattice version']                 =  PARAMS['lattice_version']
        SUMMARY['Nsigma']                          =  PARAMS['nbsigma']
        SUMMARY['injection energy [MeV]']          =  PARAMS['injection_energy']
        SUMMARY['\u03C3x_i*    [mm]']              =  PARAMS['twiss_x_i'].sigmaH()*1.e3
        SUMMARY["\u03C3x'_i* [mrad]"]              =  PARAMS['twiss_x_i'].sigmaV()*1.e3
        SUMMARY['\u03C3y_i*    [mm]']              =  PARAMS['twiss_y_i'].sigmaH()*1.e3
        SUMMARY["\u03C3y'_i* [mrad]"]              =  PARAMS['twiss_y_i'].sigmaV()*1.e3
        SUMMARY["\u03B5x_i [mrad*mm]"]             =  PARAMS['emitx_i']*1.e6
        SUMMARY["\u03B5y_i [mrad*mm]"]             =  PARAMS['emity_i']*1.e6
        SUMMARY["\u03B5w_i* [rad]"]                =  '{:8.2e} {{\u0394\u03C6,\u03B4\u03B3}}'.format(PARAMS['emitw_i'])
        SUMMARY["\u03B2w_i* [rad]"]                =  '{:8.2e} {{\u0394\u03C6,\u03B4\u03B3}}'.format(PARAMS['betaw_i'])
        SUMMARY["\u03B4\u03B3_i*"]                 =  '{:8.2e} norm. energy spread @ injection (w0)'.format(PARAMS['w0_i'])
        SUMMARY['\u0394T/T_i']                     =  '{:8.2e} kin. energy spread @ injection'.format(PARAMS['DT2T_i'])
    
    if FLAGS['dWf'] == 1:
        SUMMARY['\u0394Wmax_i* [MeV]']            =  '{:8.2e} max \u0394W on separatrix'.format(PARAMS.get('DWmax',0))
        SUMMARY['wmax*']                           =  '{:8.2e} max \u0394\u03B3 on separatrix'.format(PARAMS.get('wmax',0))
        SUMMARY['\u0394p/pmax_i [%]']             =  '{:8.2e} max \u0394p/p on separatrix'.format(PARAMS.get('Dp2pmax',0)*1.e2)
        # SUMMARY['zmax* [m]']                       =  '{:8.2e} max z on separatrix'.format(abs(PARAMS['zmax']))
        SUMMARY['\u0394p/p0*']                     =  '{:8.2e} impulse spread'.format(PARAMS['Dp2p0_i'])
        SUMMARY['\u0394\u03C60* [rad]']            =  '{:8.2e} phase'.format(PARAMS['Dphi0_i'])
        SUMMARY['|z0|* [m]']                       =  '{:8.2e} bunch length'.format(abs(PARAMS['z0_i']))
        SUMMARY['\u03C9* [MHz]']                   =  '{:8.2e} synchotron frquency'.format(PARAMS.get('omgl_0',0)*1.e-6)
    else:
        SUMMARY['separatrix:']                     =  '{}'.format('NO acceleration')
    return
def I0(x):
    """
    Modified Bessel function I of integer order 0
    ref.: Hanbook of Mathematical Functions, M.Abramowitz & I.A.Stegun
    """
    return SP.iv(0,x)
def I1(x):
    """
    Modified Bessel function I of integer order 1
    ref.: Hanbook of Mathematical Functions, M.Abramowitz & I.A.Stegun pp.379
    """
    return SP.iv(1,x)
def k0prot(gradient=0.,tkin=0.):         #TODO still used?
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
def scalek0prot(k0=0.,tki=0.,tkf=0.):    #TODO still used?
    """Scale Quadrupole strength k0 for increase of kin. energy from tki to tkf  (only for protons!)"""
    bgi  = Proton(tki).gamma_beta
    bgf  = Proton(tkf).gamma_beta
    k= k0 * bgi/bgf
    return k
def dBdxprot(k0=0.,tkin=0.):             #TODO still used?
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
    # print dictionary sorted by key
    print(dictstring(what,text,filter,njust),end=end)
def print_verbose(level,*args):
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
class Test_set_utilities(unittest.TestCase):
    def test_particle(self):
        print('----------------------------------test_particle')
        print('Proton(5.)\n'+Proton(5.).toString())
        print('Electron(5.)\n'+Electron(5.).toString())
    def test_params(self):
        print('----------------------------------test_params')
        dictprnt( PARAMS,text=' PARAMS')
    def test_particle_energy_adjustment(self):
        print('----------------------------------test_particle_energy_adjustment')
        p = Proton(50.)
        print(repr(p)+':')
        print(p.toString())
        p1 = p(100.)
        print(repr(p1)+':')
        print(p1.toString())
        try:
            p2 = p(-10.)
        except ValueError:
            print("This was a test: Particle energy deliberately set negative!")
    def test_ellipse(self):
        print('----------------------------------test_ellipse')
        def ellicp(xy,alfa,beta,emit):
            """ convert twiss parameters to plot parameters """
            gamma = (1.+alfa**2)/beta
            H = 0.5*(beta+gamma)     # see CERN's Formelsammlung
            a = sqrt(0.5*emit)*(sqrt(H+1.)+sqrt(H-1.))
            b = sqrt(0.5*emit)*(sqrt(H+1.)-sqrt(H-1.))
            tilt = degrees(0.5*atan(2*alfa/(gamma-beta)))
            # return plot prameters as  (origin,width,height,tilt)
            return (xy,a,b,tilt)
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
    def test_lissajous_plot(self):
        print('----------------------------------test_lissajous_plot')
        names = ('t','SIN','COS')
        fn = Functions(names)
        for ph in range(0,361):
            phr = radians(ph)
            fn.append(ph,(sin(phr),cos(3*phr)))
        
        t= [fn(i,'t')   for i in range(fn.nbpoints)]
        S= [fn(i,'SIN') for i in range(fn.nbpoints)]
        C= [fn(i,'COS') for i in range(fn.nbpoints)]
        
        for i in range(fn.nbpoints): 
            if 0: print(t[i],S[i],C[i])
        plt.plot(S,C)
        plt.show()
if __name__ == '__main__':
    unittest.main()
