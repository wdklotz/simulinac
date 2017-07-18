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
from math import pi,sqrt,sin, cos, radians, degrees
import logging

## create logger
logger = logging.getLogger("logger")
logger.setLevel(logging.DEBUG)
#create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
#create formatter
formatter = logging.Formatter("%(levelname)s: %(filename)s[%(lineno)d] %(message)s")
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)
## DEBUG
def DEBUG(string,arg=''):
    """
    Print debug message
    IN:
        string to print 
        values to print
    """
    if isinstance(arg,list):
        print('DEBUG: {} \nlist={}'.format(string,arg))
    elif isinstance(arg,dict):
        print('DEBUG: {} \ndict={}'.format(string,arg))
    else:
        print('DEBUG: ',string,arg)
## defaults
class Defaults(object):
    def __init__(self):
        self.priv = {
            'lichtgeschwindigkeit': 299792458.,    # [m/s] const
            'elementarladung': 1.602176565e-19,    # [coulomb] const
            'proton_mass': 938.272,      # [MeV/c**2] const
            'electron_mass': 0.5109989,  # [MeV/c**2] const
            }
        self.conf = {            ## CONFIG constants and setutil ...
            'Ez_feld': 1.04,             # [MV/m] default
            'spalt_laenge': 0.02,        # [m] default
            'cavity_laenge': 0.08,       # [m] default
            'soll_phase': -10.0,         # [deg] default
            'frequenz': 813.e6,          # [Hz] default
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

CONF = Defaults()
CONF['wellenlänge']     = CONF['lichtgeschwindigkeit']/CONF['frequenz']
CONF['Dz']              = CONF['wellenlänge']/18.  # Dz aka delta-z is 1/18-th of wavelength per default
CONF['spalt_spannung']  = CONF['Ez_feld']*CONF['spalt_laenge']

## relativistic particle
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
## proton
class Proton(Particle):
    def __init__(self,tkin=CONF['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass=CONF['proton_mass'],name='proton')
## electron
class Electron(Particle):
    def __init__(self,tkin=CONF['injection_energy']):
        super(Electron,self).__init__(tkin=tkin,mass=CONF['electron_mass'],name='electron')
## the default reference particle
# Particle.soll = Proton()
CONF['sollteilchen'] = Proton()
## utilities
def epsiz(dz,beta,gamma,trtf):
    """
    Helper to calculate longitudinal phase space ellipse parameters
    Ellipse nach T.Wangler (6.47) pp.185
    (w/w0)**2 + (dphi/dhi0)**2 = 1 entspricht (gamma*x)**2 + (beta*x')**2 = epsilon
    """
    qE0   = CONF['Ez_feld']
    lamb  = CONF['wellenlänge']
    m0c2  = CONF['proton_mass']
    dphi0 = - radians((360./beta/lamb)*dz)     #umrechnung [m] -> [rad]
    R = qE0*lamb*sin(-radians(CONF['soll_phase']))
    R = R * beta*beta*beta*gamma*gamma*gamma*trtf
    R = R / (2.*pi*m0c2)
    R = sqrt(R)
    w0     = R * dphi0
    epsi   = w0 * dphi0
    sigw   = sqrt(epsi*R)
    sigphi = -sqrt(epsi/R)    #[rad]
    sigDp  = gamma/(gamma+1.)*sigw      #umrechnung DW/W -> Dp/p, W kin. energie
    sigDz  = - beta*lamb/360.*degrees(sigphi)   #[m]
    sigphi = degrees(sigphi)
    return dict(epsi=epsi,dz=dz,sigDz=sigDz,sigDp=sigDp,sigw=sigw,sigphi=sigphi)

result = epsiz(CONF['Dz'],
                CONF['sollteilchen'].beta,
                CONF['sollteilchen'].gamma,
                CONF['sollteilchen'].trtf(CONF['spalt_laenge'],CONF['frequenz'])
                )

CONF['emitz_i']   = result['epsi']
CONF['Dp/p']      = result['sigDp']
CONF['sigPhiz_i'] = result['sigphi']
CONF['DW/W']      = result['sigw']

SUMMARY = {}

def collect_summary():
    SUMMARY['frequency [Hz]'] = CONF['frequenz']
    SUMMARY['QF gradient [T/m]'] = CONF['qf_gradient']
    SUMMARY['QD gradient [T/m]'] = CONF['qd_gradient']
    SUMMARY['Quad pole length [m]'] = CONF['ql']
    SUMMARY['Quad bore radius [m]'] = CONF['quad_bore_radius']
    SUMMARY['injection energy [MeV]'] = CONF['injection_energy']
    SUMMARY['emitx_i [rad*m]'] = CONF['emitx_i']
    SUMMARY['emity_i [rad*m]'] = CONF['emity_i']
    SUMMARY['emitz_i* [rad]'] = CONF['emitz_i']
    SUMMARY['sigx_i* [mm]'] = 1000.*sqrt(CONF['betax_i']*CONF['emitx_i'])  # enveloppe @ entrance
    SUMMARY['sigy_i* [mm]'] = 1000.*sqrt(CONF['betay_i']*CONF['emity_i'])
    SUMMARY['sync. phase [deg]'] = CONF['soll_phase']
    SUMMARY['cavity gap length [m]'] = CONF['spalt_laenge']
    SUMMARY['cavity length [m]'] = CONF['cavity_laenge']
    SUMMARY['wavelength [m]'] = CONF['wellenlänge']
    SUMMARY['cavity gap voltage* [MV]'] = CONF['spalt_spannung']
    SUMMARY['acc. field Ez [MV/m]'] = CONF['Ez_feld']
    SUMMARY['lattice version'] = CONF['lattice_version']
    SUMMARY['QF pole strength* [T]'] = CONF['qf_gradient'] * CONF['ql']
    SUMMARY['QF current* [A/winding]'] = (CONF['qf_gradient'] * (CONF['ql']*1000.)**2 )/2.52/CONF['n_coil']
    SUMMARY['QF power estimate* [W]'] = 0.0115 *SUMMARY['QF current* [A/winding]']**2  # R=0.0115 Ohms
    SUMMARY['QF coil [windings]'] = CONF['n_coil']
    SUMMARY['Dz_i(bunch spread) [m]'] = CONF['Dz']
    SUMMARY['Dp/p_i(impulse spread)* [%]'] = CONF['Dp/p']*1.e+2
    SUMMARY['Dph/ph_i(phase spread)* [deg]'] = CONF['sigPhiz_i']
    SUMMARY['DW/m0c2_i(energy spread)* [%]'] = CONF['DW/W']*1.e+2
    SUMMARY['DW/m0c2 max* [%]'] = wakzp = accpt_w(    # energy acceptance in %
            CONF['Ez_feld'],
            CONF['sollteilchen'].trtf(CONF['spalt_laenge'],CONF['frequenz']),
            CONF['soll_phase'],
            CONF['wellenlänge'],
            CONF['sollteilchen'])*1.e+2
    SUMMARY['Dp/p max* [%]'] = 1./(1.+1./CONF['sollteilchen'].gamma)*wakzp  # impule acceptanc in %
    return

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
        DEBUG('accpt_w: tk {:4.4f} trtf {:4.4f} res {:4.4e}'.format(tk,trtf,res))
        phsoll = radians(phis)
        res = res*(phsoll*cos(phsoll)-sin(phsoll))
        DEBUG('(dW/W)^2 nach T.Wangler pp. 179: res',res)
        res = sqrt(res)/(particle.gamma-1.)
    else:
        res = 0.
    return res

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
    Helper to print objects as dictionary
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

def dictprnt(what,text='',filter=[],njust=30): 
    """
    Helper to print dictionaries
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
    for k,v in sorted(what.items()):
        if k in filter:
            continue
        print(k.rjust(njust),':',v)
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
    Print a nice table to memory.
    IN:
        headr = table header [...]
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
    result = epsiz(CONF['Dz'],
              CONF['sollteilchen'].beta,
              CONF['sollteilchen'].gamma,
              CONF['sollteilchen'].trtf(CONF['spalt_laenge'],CONF['frequenz'])
              )
    for k,v in result.items():
        print('{}\t{:g}'.format(k,v))
## main
if __name__ == '__main__':
    test0()
    # test1()

