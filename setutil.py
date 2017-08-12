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
from math import pi,sqrt,sin, cos, radians, degrees, pow
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
        pp   = pprint.PrettyPrinter(indent=4)  ## use pprint module
        sarg = pp.pformat(arg)
        print('DEBUG: {} typ(dict)\n{}'.format(string,sarg))
    else:
        print('DEBUG: {}{}'.format(string,arg))

class Defaults(object):
    def __init__(self):
        self.priv = {
            'lichtgeschwindigkeit': 299792458.,    # [m/s] const
            'elementarladung': 1.602176565e-19,    # [coulomb] const
            'proton_mass': 938.272,      # [MeV/c**2] const
            'electron_mass': 0.5109989,  # [MeV/c**2] const
            }
        self.conf = {            ## CONFIG constants and setutil ...
            'Ez_feld': 1.00,             # [MV/m] default
            'spalt_laenge': 0.02,        # [m] default
            'cavity_laenge': 0.08,       # [m] default
            'soll_phase': -30.0,         # [deg] default
            'frequenz': 814.e6,          # [Hz] default
            'injection_energy': 50.,     # [MeV] default
            'qf_gradient': 16.0,         # [T/m] default
            'qd_gradient': 16.0,         # [T/m] default
            'quad_bore_radius': 0.02,    # Vorgabe quadrupole bore radius [m]
            'emitx_i' : 1.e-6,           # [m*rad] Vorgabe emittance @ entrance
            'emity_i' : 1.e-6,           # [m*rad] Vorgabe emittance @ entrance
            'emitz_i' : 7.7e-4,          # longitudinal emittance T.Wangler (6.49) pp.186
            'betax_i' : 0.780,           # [m] Vorgabe twiss betax @ entrance
            'betay_i' : 2.373,           # [m] Vorgabe twiss betax @ entrance
            'alfax_i' : 0.0,             # Vorgabe twiss alphax @ entrance
            'alfay_i' : 0.0,             # Vorgabe twiss alphaxy @ entrance
            'Dz'      : 0.02,            # [m] Vorgabe longitudinal displacement Dz
            'Dp/p'    : 3.55654e-4,      # [rad] Vorgabe relative impulse Dp/p
            'dWf'     : False,           # acceleration on/off flag default
            'periodic': True,            # periodic lattice? default
            'verbose' : 1,               # print flag (True) default
            'egf'     : False,           # emittance grow flag default
            'sigma'   : True,            # beam sizes by sigma-tracking
            'n_coil'  : 30               # nbof coil windings
            }
    def __getitem__(self,key):
        if key in self.conf:
            return self.conf[key]
        elif key in self.priv:
            return self.priv[key]
        else:
            return None
    def __setitem__(self,key,value):
        self.conf[key]=value
    def items(self):
        res={}
        for k,v in self.priv.items():
            res[k]=v
        for k,v in self.conf.items():
            res[k]=v
        return res.items()

## BLOCKDATA "CONF"
CONF = Defaults()       # global data block called CONF (like Fortran's BLOCKDATA)
CONF['wellenlänge']     = CONF['lichtgeschwindigkeit']/CONF['frequenz']
CONF['Dz']              = CONF['wellenlänge']/36.  # Dz (aka delta-z) is 1/36-th of wavelength (i.e.10 deg per default)
CONF['spalt_spannung']  = CONF['Ez_feld']*CONF['spalt_laenge']

class Particle(object):                          
    # soll = None  # class member: reference particle a.k.a. soll Teilchen - deactivated, caused serious error
    def __init__(self,tkin=0.,mass=CONF['proton_mass'],name='proton'):
        self._set_self(tkin,mass,name)
    def _set_self(self,tkin,mass,name):
        self.tkin       = tkin                     # kinetic energy [MeV]
        self.e0         = mass                     # rest mass [MeV/c**2]
        self.e          = self.e0+self.tkin        # total energy [MeV]
        self.gamma      = self.e/self.e0
        self.beta       = sqrt(1.-1./(self.gamma*self.gamma))
        self.gamma_beta = self.gamma * self.beta
        self.p          = self.gamma_beta * self.e0                 # impulse [Mev/c]
        self.v          = self.beta*CONF['lichtgeschwindigkeit']    # velocity [m/s]
        self.brho       = 1.e+6/CONF['lichtgeschwindigkeit']*self.gamma_beta*self.e0 # [T*m]
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
        teta = 2.*pi*fRF*gap / (self.beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf
    def __call__(self,tkin):  # call Particle instance to change its kin. energy
        self._set_self(tkin=tkin,mass=self.e0,name=self.name)
        return self

class Proton(Particle):
    def __init__(self,tkin=CONF['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass=CONF['proton_mass'],name='proton')

class Electron(Particle):
    def __init__(self,tkin=CONF['injection_energy']):
        super(Electron,self).__init__(tkin=tkin,mass=CONF['electron_mass'],name='electron')

## Sollteichen
CONF['sollteilchen'] = Proton()

## long. emittance
def epsiz(particle=CONF['sollteilchen'],gap=0.0,trtf=0.75):
    """
    Helper to calculate longitudinal phase space ellipse parameters
    Ellipse nach T.Wangler (6.47) pp.185
    (w/w0)**2 + (dphi/dhi0)**2 = 1 entspricht (gamma*x)**2 + (beta*x')**2 = epsilon, alpha=0
    """
    # NOTE: epsiz should be an atrtribute of RF cavities (todo?)
    qE0      = CONF['Ez_feld']
    lamb     = CONF['wellenlänge']
    m0c2     = particle.e0
    dz       = CONF['Dz']
    gb       = particle.gamma_beta
    beta     = particle.beta
    gamma    = particle.gamma
    if gap != 0.: trtf = particle.trtf(gap,CONF['frequenz'])
    dphi0    = -radians((360./beta/lamb)*dz)     # Umrechnung [m] -> [rad]
    
    R = qE0*lamb*sin(-radians(CONF['soll_phase']))  ## R ist die grosse Wurzel
    R = R*pow(gb,3)*trtf
    R = R/(2.*pi*m0c2)
    R = sqrt(R)
    
    w0     = R*dphi0
    epsi   = w0*dphi0
    sigw   = sqrt(epsi*R)
    sigphi = -sqrt(epsi/R)    #[rad]
    sigphi = degrees(sigphi)  #[deg]
    sigDp  = gamma/(gamma+1.)*sigw      # Umrechnung DW/W -> Dp/p, W kin. Energie
    sigDz  = -beta*lamb/360.*sigphi    # [m]

    CONF['emitz_i']   = epsi
    CONF['Dp/p']      = sigDp
    CONF['sigphi']    = sigphi
    CONF['DW/W']      = sigw
    return dict(epsi=epsi,dz=dz,sigDz=sigDz,sigDp=sigDp,sigw=sigw,sigphi=sigphi)
epsiz()     ## calculate the long. emittance with def. parameters

def accpt_w(Ez,trtf,phis,lamb,particle):
    """
    Energieakzeptanz dW/W nach T.Wangler pp. 179
    Ez       - [Mev/m] accelerating field gradient on gap axis
    trtf     - transit time factor
    phis     - [deg] sync. phase (-90 to 0)
    lamb     - [m] rf wave length
    particle - sync. particle
    """
    if CONF['dWf']:
        gb = particle.gamma_beta
        m0 = particle.e0  # rest mass
        tk = particle.tkin
        res = 2.*Ez*trtf*gb*gb*gb*lamb/pi/m0
        # DEBUG('accpt_w: tk {:4.4f} trtf {:4.4f} res {:4.4e}'.format(tk,trtf,res))
        phsoll = radians(phis)
        res = res*(phsoll*cos(phsoll)-sin(phsoll))
        # DEBUG('(dW/W)^2 nach T.Wangler pp. 179: res',res)
        res = sqrt(res)/(particle.gamma-1.)
    else:
        res = 0.
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

    sections = CONF['sections']   ## comes from INPUT
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

    SUMMARY['sections']                        = CONF['sections']
    SUMMARY['frequency [Hz]']                  = CONF['frequenz']
    SUMMARY['Quad bore radius [m]']            = CONF['quad_bore_radius']
    SUMMARY['injection energy [MeV]']          = CONF['injection_energy']
    SUMMARY['(emitx)i [rad*m]']                = CONF['emitx_i']
    SUMMARY['(emity)i [rad*m]']                = CONF['emity_i']
    SUMMARY['(emitz)i* [rad]']                 = CONF['emitz_i']
    SUMMARY['(sigx)i* [mm]']                   = 1000.*sqrt(CONF['betax_i']*CONF['emitx_i'])  # enveloppe @ entrance
    SUMMARY['(sigy)i* [mm]']                   = 1000.*sqrt(CONF['betay_i']*CONF['emity_i'])
    SUMMARY['wavelength [m]']                  = CONF['wellenlänge']
    SUMMARY['lattice version']                 = CONF['lattice_version']
    SUMMARY['(Dz)i [m]']                       = CONF['Dz']
    SUMMARY['(Dp/p)i(impulse spread)* [%]']    = CONF['Dp/p']*1.e+2
    SUMMARY['(Dphi)i* [deg]']                  = CONF['sigphi']
    SUMMARY['(DW/m0c2)i(energy spread)* [%]']  = CONF['DW/W']*1.e+2
    SUMMARY['(DW/m0c2)accept limit* [%]']      = wakzp = accpt_w(CONF['Ez_feld'],CONF['sollteilchen'].trtf(CONF['spalt_laenge'],CONF['frequenz']),CONF['soll_phase'],CONF['wellenlänge'],CONF['sollteilchen'])*1.e+2  ## energy acceptance
    SUMMARY['(Dp/p)accept limit* [%]']         = 1./(1.+1./CONF['sollteilchen'].gamma)*wakzp  # impule acceptanc in %
    return

## utilities
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
    verbose = CONF['verbose']
    if verbose >= level:
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
    dictprnt(CONF,text='CONF')    
    
    print('Sollteilchen\n'+CONF['sollteilchen'].string())
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
    print('test epsiz(): the helper to calculate longitudinal phase space parameters')
    result = epsiz()
    for k,v in result.items():
        print('{}\t{:g}'.format(k,v))

## main
if __name__ == '__main__':
    test0()
    test1()

