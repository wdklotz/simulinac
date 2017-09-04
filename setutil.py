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
from math import pi,sqrt,sin,cos,radians,degrees,pow,fabs,exp
import logging, pprint

## logger
ch        = logging.StreamHandler()     ## console handler
ch.setLevel(logging.DEBUG)              ## set handler level
formatter = \
    logging.Formatter("%(levelname)s: %(filename)s[%(lineno)d] %(message)s")
ch.setFormatter(formatter)              ## set handler's format
logger    = logging.getLogger("logger")
logger.addHandler(ch)                   ## add handler to logger

def DEBUG(string,arg=''):
    """
    Print debug message
    IN:
        string: text to print 
        arg:  values to print
    """
    if isinstance(arg,list):
        # print('DEBUG: {} \nlist={}'.format(string,arg))
        pp   = pprint.PrettyPrinter(indent=4)  ## use pprint module
        sarg = pp.pformat(arg)
        print('DEBUG: {} typ(list)\n{}'.format(string,sarg))
    elif isinstance(arg,dict):
        # print('DEBUG: {} \ndict={}'.format(string,arg))
        pp   = pprint.PrettyPrinter(indent=4,width=60)  ## use pprint module
        sarg = pp.pformat(arg)
        print('DEBUG: {} typ(dict)\n{}'.format(string,sarg))
    else:
        print('DEBUG: {}{}'.format(string,arg))


## DEFAULTS "FLAGS" & "PARAMS"
FLAGS  = dict(
        periodic             = False,            # periodic lattice? default
        egf                  = True,             # emittance grow flag default
        sigma                = True,             # beam sizes by sigma-tracking
        KVprint              = False,            # print a dictionary of Key-Value pairs, no display
        map                  = True,             # use maps to track trajectories through RFGap
        dWf                  = 1.,               # acceleration on/off flag 1=on,0=off
        verbose              = 0                 # print flag default = 0
        )
PARAMS = dict(
        lichtgeschwindigkeit = 299792458.,       # [m/s] const
        elementarladung      = 1.602176565e-19,  # [coulomb] const
        proton_mass          = 938.272,          # [MeV/c**2] const
        electron_mass        = 0.5109989,        # [MeV/c**2] const
        Ez_feld              = 1.00,             # [MV/m] default
        spalt_laenge         = 0.02,             # [m] default
        cavity_laenge        = 0.08,             # [m] default
        soll_phase           = -30.0,            # [deg] default
        frequenz             = 814.e6,           # [Hz] default
        injection_energy     = 50.,              # [MeV] default
        qf_gradient          = 16.0,             # [T/m] default
        qd_gradient          = 16.0,             # [T/m] default
        quad_bore_radius     = 0.02,             # Vorgabe quadrupole bore radius [m]
        n_coil               = 30,               # nbof coil windings
        emitx_i              = 1.e-6,            # [m*rad] Vorgabe emittance @ entrance
        emity_i              = 1.e-6,            # [m*rad] Vorgabe emittance @ entrance
        betax_i              = 0.780,            # [m] Vorgabe twiss betax @ entrance
        betay_i              = 2.373,            # [m] Vorgabe twiss betax @ entrance
        alfax_i              = 0.0,              # Vorgabe twiss alphax @ entrance
        alfay_i              = 0.0,              # Vorgabe twiss alphaxy @ entrance
        sigmaz_i             = 0.02,             # [m] max long. half-width displacement
        aperture             = 0.011,            # aperture = bore radius
        )
PARAMS['wellenlänge']     = PARAMS['lichtgeschwindigkeit']/PARAMS['frequenz']
PARAMS['sigmaz_i']        = PARAMS['wellenlänge']/36.  # sigma-z is 1/36-th of wavelength (i.e.10 deg per default)
PARAMS['spalt_spannung']  = PARAMS['Ez_feld']*PARAMS['spalt_laenge']

CpValues = dict(z=0.,sigma_x=0.,sigma_y=0.,Tkin=0.)

class Particle(object):                          
    # soll = None  # class member: reference particle a.k.a. soll Teilchen - deactivated, caused serious error
    def __init__(self,tkin=0.,mass= PARAMS['proton_mass'],name='proton'):
        self._set_self(tkin,mass,name)
    def _set_self(self,tkin,mass,name):
        self.tkin       = tkin                     # kinetic energy [MeV]
        self.e0         = mass                     # rest mass [MeV/c**2]
        self.e          = self.e0+self.tkin        # total energy [MeV]
        self.gamma      = self.e/self.e0
        self.beta       = sqrt(1.-1./(self.gamma*self.gamma))
        self.gamma_beta = self.gamma * self.beta
        self.p          = self.gamma_beta * self.e0                 # impulse [Mev/c]
        self.v          = self.beta* PARAMS['lichtgeschwindigkeit']    # velocity [m/s]
        self.brho       = 1.e+6/ PARAMS['lichtgeschwindigkeit']*self.gamma_beta*self.e0 # [T*m]
        self.name       = name
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
    def __call__(self,tkin):  # call Particle instance to change its kin. energy
        self._set_self(tkin=tkin,mass=self.e0,name=self.name)
        return self

class Proton(Particle):
    def __init__(self,tkin= PARAMS['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass= PARAMS['proton_mass'],name='proton')

class Electron(Particle):
    def __init__(self,tkin= PARAMS['injection_energy']):
        super(Electron,self).__init__(tkin=tkin,mass= PARAMS['electron_mass'],name='electron')

## Sollteichen
PARAMS['sollteilchen'] = Proton()

## long. emittance
def zellipse(sigmaz,qE0,lamb,phis,gap,particle):
    """
    Helper to calculate longitudinal phase space ellipse parameters
    Ellipse nach T.Wangler (6.47) pp.185
    (w/w0)**2 + (dphi/dhi0)**2 = 1 entspricht gammaz*w0**2 + betaz*Dphi0**2 = emittz, alphaz=0
    IN:
        sigmaz:   max half-width long. displacement [m]
        qE0:      gap field strength [MV/m]
        lamb:     wavelength [m]
        phis:     syncronous phase [rad]
        gap:      cavity gap [m]
        particle: Particle object
    """
    # NOTE: zellipse should be an atrtribute of RF cavities (I am not sure!?)
    m0c2      = particle.e0
    gb        = particle.gamma_beta
    beta      = particle.beta
    gamma     = particle.gamma
    if gap != 0.: 
        T = particle.trtf(gap, PARAMS['frequenz'])
    else:
        T = 0.75       # have to take some reasonable value!
    
    # large amplitude oscillations (T.Wangler pp. 175)
    wmax  = sqrt(2.*qE0*T*pow(gb,3)*lamb/(pi*m0c2)*(phis*cos(phis)-sin(phis)))  # T.Wangler (6.28)
    DWmax = wmax*m0c2       # [MeV] conversion ---> [MeV]
    phi1s = -phis           # [rad]
    phi2s = 2.*phis         # [rad] Naehrung T.Wangler pp.178
    psis  = 3.*fabs(phis)   # [rad] Naehrung T.Wangler pp.178
    
    # small amplitude oscillations (T.Wangler pp.184)
    kl02 = 2.*pi*qE0*T*sin(-phis)/(m0c2*pow(gb,3)*lamb)
    omegal0 = sqrt(kl02)*beta* PARAMS['lichtgeschwindigkeit']
    omegal0_div_omega = sqrt(qE0*T*lamb*sin(-phis)/(2.*pi*m0c2*pow(gamma,3)*beta))
    
    # Dphi0 = (phi0 - phis) maximum half-width phase dispersion see T.Wangler
    Dphi0  = (2.*pi*sigmaz)/(beta*lamb)     # [rad]  conv. z --> phi
    w0     = sqrt(qE0*T*pow(gb,3)*lamb*sin(-phis)*pow(Dphi0,2)/(2.*pi*m0c2))
    DW     = w0*m0c2              # conversion --> [MeV]
    emitz  = Dphi0*DW             # [rad*MeV]
    gammaz = pow(Dphi0,2)/emitz   # [rad/MeV]
    betaz  = pow(DW,2)/emitz      # [MeV/rad]
    alphaz = 0.                   # always!

    res =  dict(
            Dphi0           = Dphi0,
            DWmax           = DWmax,
            Dphimax         = psis,
            omegal0         = omegal0,
            w0              = w0,
            emitz           = emitz,
            gammaz          = gammaz,
            betaz           = betaz,
            DW              = DW,
            alphaz          = alphaz )
    res['omegal0/omega']    = omegal0_div_omega
    res['Dphimax+']         = phi1s
    res['Dphimax-']         = phi2s
    return res

## data for summary
SUMMARY = {}
def collect_data_for_summary(lattice):
    def elements_in_lattice(typ,sec):
        '''
        Filter elements of class <typ> and section <sec> form lattice
        IN:
            lattice = object [Lattice]
            typ     = element class [string]
            sec     = section name [string]
        OUT:
            iterator of filtered elements
        '''
        import itertools
        def predicate(element):
            try:
                test = (type(element[0]).__name__ == typ and element[0].sec == sec)
            except AttributeError:
                test = (type(element[0]).__name__ == typ)  ## no sec tags? take all!
            return not test
        filtered_elements = itertools.filterfalse(predicate,lattice.seq)
        # for t in filtered_elements: DEBUG('filterfalse',(t,t[0].label))  ## whazit
        return filtered_elements

    def elements_in_section(typ,sec):
        """
        Get a list of elements of same type in a section
        """
        elements = list(elements_in_lattice(typ,sec))
        new_elements = []
        seen = set()                  ## helper to eliminate duplicate entries
        for itm in elements:
            label = itm[0].label      ## itm is tupel (element,s0,s1)
            if label in seen:
                continue
            else:
                seen.add(label)
                new_elements.append(itm[0])
        return new_elements

    sections =  PARAMS['sections']   ## comes from INPUT
    if len(sections) == 0: sections = ['*']      ## section wildcart
    types = ['QF','QD']
    for sec in sections:
        for typ in types:
            elements = elements_in_section(typ,sec)
            for itm in elements:
                k0 = itm.k0
                dBdz = k0*itm.particle.brho
                length = itm.length
                # SUMMARY['{2} [{1}.{0}]    k0 [m^-2]'.format(sec,typ,itm.label)] = k0
                SUMMARY['{2} [{1}.{0}]   dBdz [T/m]'.format(sec,typ,itm.label)] = dBdz
                SUMMARY['{2} [{1}.{0}]   length [m]'.format(sec,typ,itm.label)] = length

                PARAMS['{2}[{1}.{0}]dBdZ'.format(sec,typ,itm.label)] = dBdz
                PARAMS['{2}[{1}.{0}]length'.format(sec,typ,itm.label)] = length
    types = ['RFG']
    for sec in sections:
        for typ in types:
            elements = elements_in_section(typ,sec)
            for itm in elements:
                gap     = itm.gap
                Ez      = itm.u0/gap
                PhiSoll = degrees(itm.phis)
                # length  = itm.length
                SUMMARY['{2} [{1}.{0}]  gap    [m]'.format(sec,typ,itm.label)] = gap
                SUMMARY['{2} [{1}.{0}]  Ez  [MV/m]'.format(sec,typ,itm.label)] = Ez
                SUMMARY['{2} [{1}.{0}]  phis [deg]'.format(sec,typ,itm.label)] = PhiSoll

                PARAMS['{2}[{1}.{0}]gap'.format(sec,typ,itm.label)] = gap
                PARAMS['{2}[{1}.{0}]Ez'.format(sec,typ,itm.label)] = Ez
                PARAMS['{2}[{1}.{0}]phis'.format(sec,typ,itm.label)] = PhiSoll
    types = ['RFC']
    for sec in sections:
        for typ in types:
            elements = elements_in_section(typ,sec)
            for itm in elements:
                gap     = itm.gap
                Ez      = itm.u0/gap
                PhiSoll = degrees(itm.phis)
                length  = itm.length
                SUMMARY['{2} [{1}.{0}]  gap    [m]'.format(sec,typ,itm.label)] = gap
                SUMMARY['{2} [{1}.{0}]  Ez  [MV/m]'.format(sec,typ,itm.label)] = Ez
                SUMMARY['{2} [{1}.{0}]  phis [deg]'.format(sec,typ,itm.label)] = PhiSoll
                SUMMARY['{2} [{1}.{0}]  length [m]'.format(sec,typ,itm.label)] = length

                PARAMS['{2}[{1}.{0}]gap'.format(sec,typ,itm.label)] = gap
                PARAMS['{2}[{1}.{0}]Ez'.format(sec,typ,itm.label)] = Ez
                PARAMS['{2}[{1}.{0}]phis'.format(sec,typ,itm.label)] = PhiSoll
                PARAMS['{2}[{1}.{0}]length'.format(sec,typ,itm.label)] = length

    SUMMARY['track with map']                  =  FLAGS['map']
    SUMMARY['sigma tracking']                  =  FLAGS['sigma']
    SUMMARY['emittance growth']                =  FLAGS['egf']
    SUMMARY['ring lattice']                    =  FLAGS['periodic']
    SUMMARY['accON']                           =  False if  FLAGS['dWf'] == 0. else  True
    SUMMARY['frequency [MHz]']                 =  PARAMS['frequenz']*1.e-6
    SUMMARY['Quad bore radius [m]']            =  PARAMS['quad_bore_radius']
    SUMMARY['injection energy [MeV]']          =  PARAMS['injection_energy']
    SUMMARY['(emitx)i [mrad*mm]']              =  PARAMS['emitx_i']*1.e6
    SUMMARY['(emity)i [mrad*mm]']              =  PARAMS['emity_i']*1.e6
    SUMMARY['(emitz)i* [rad*KeV]']             =  PARAMS['emitz_i']*1.e3
    SUMMARY['(sigx)i* [mm]']                   =  sqrt( PARAMS['betax_i']* PARAMS['emitx_i'])*1.e3  # enveloppe @ entrance
    SUMMARY['(sigy)i* [mm]']                   =  sqrt( PARAMS['betay_i']* PARAMS['emity_i'])*1.e3
    SUMMARY['wavelength* [cm]']                =  PARAMS['wellenlänge']*1.e2
    SUMMARY['lattice version']                 =  PARAMS['lattice_version']
    SUMMARY['(sigmaz)i [mm]']                  =  PARAMS['sigmaz_i']*1.e3
    SUMMARY['(DW)i* [KeV]']                    =  PARAMS['DW']*1.e3
    SUMMARY['(Dphi)i* [deg]']                  =  degrees( PARAMS['Dphi0'])
    SUMMARY['(DW)max* [KeV]']                  =  PARAMS['DWmax']*1.e3        # energy acceptance
    SUMMARY['max bunch length* [deg]']         =  degrees( PARAMS['Dphimax'])  # phase acceptance
    SUMMARY['betaz_i* [KeV/rad]']              =  PARAMS['betaz_i']*1.e3
    SUMMARY['gammaz_i* [rad/KeV]']             =  PARAMS['gammaz_i']*1.e-3
    SUMMARY['synchrotron freq_i* [MHz]']       =  PARAMS['omegal0']*1.e-6
    SUMMARY['sync.freq_i/rf_freq* [%]']        =  PARAMS['omegal0/omega']*1.e2
    SUMMARY['aperture [m]']                    =  PARAMS['aperture']
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
        # DEBUG('I0->x ',x)
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
        # DEBUG('I0->x ',x)
        res = res*exp(x)/sqrt(x)
    else:
        raise RuntimeError('I0(): x={} negative argument!'.format(x))
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
        res = res*exp(x)/sqrt(x)
    else:
        raise RuntimeError('I1(x): negative argument!')
        sys.exit(1)
    return res
## utilities
def do_actions(actions):
    if 'sigma_x' in actions:
        # print('(sigma)x @ z {:8.4f}[m] = {:8.4f}[mm]'.format(CpValues['z'],CpValues['sigma_x']*1.e3))
        SUMMARY['sigma-x({:8.4f}[m])[mm]'.format(CpValues['z'])] = CpValues['sigma_x']*1.e3
        PARAMS['sigma-x({:0=6.2f})'.format(CpValues['z'])]=CpValues['sigma_x']*1.e3
    if 'sigma_y' in actions:
        # print('(sigma)y @ z {:8.4f}[m] = {:8.4f}[mm]'.format(CpValues['z'],CpValues['sigma_y']*1.e3))
        SUMMARY['sigma-y({:8.4f}[m])[mm]'.format(CpValues['z'])] = CpValues['sigma_y']*1.e3
        PARAMS['sigma-y({:0=6.2f})'.format(CpValues['z'])]=CpValues['sigma_y']*1.e3
    # if 'Tkin' in actions:
    #     print('Tkin')

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
        raise RuntimeError('setutil.k0(): negative kinetic energy?')
        sys.exit(1)

def scalek0prot(k0=0.,tki=0.,tkf=0.):
    """
    scale Quadrupole strength k0 for increase of kin. energy from tki to tkf  (only for protons!)
    """
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
        raise RuntimeError('setutil.k0(): negative kinetic energy?')
        sys.exit(1)

def objprnt (what,text='',filter=[]):
    """
    Custom helper to print objects as dictionary
    """
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
        print(k.rjust(30),':',v)
    return

def dictprnt(what,text='',filter=[],njust=35): 
    """
    Custom helper to print dictionaries (clever!?)
    """
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
    """
    Multilevel printing using verbose flag
    """
    verbose = FLAGS['verbose']
    if verbose >= level and not FLAGS['KVprint']:
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
    """
    Simple matrix print
    """
    s = [['{:+.3e}  '.format(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = ''.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    return '\n'.join(table)

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
    print('test zellipse(): the helper to calculate longitudinal phase space parameters')
    result = zellipse(
            PARAMS['sigmaz_i'],
            PARAMS['Ez_feld'],
            PARAMS['wellenlänge'],
            radians(PARAMS['soll_phase']),
            PARAMS['spalt_laenge'],
            PARAMS['sollteilchen']
            )
    for k,v in result.items():
        print('{}\t{:g}'.format(k,v))

## main
if __name__ == '__main__':
    test0()
    test1()

