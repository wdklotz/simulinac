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

MDIM=10

CONF = {                                      ## CONFIG constants and setup ...
    'lichtgeschwindigkeit': 299792458.,    # [m/s] const
    'elementarladung': 1.602176565e-19,    # [coulomb] const
    'proton_mass': 938.272,      # [MeV/c**2] const
    'electron_mass': 0.5109989,  # [MeV/c**2] const
    'Ez_feld': 1.04,             # [MV/m] default
    'spalt_laenge': 0.02,        # [m] default
    'cavity_laenge': 0.08,       # [m] default
    'soll_phase': -10.0,         # [deg] default
    'frequenz': 813.e6,          # [Hz] default
    'injection_energy': 50.,     # [MeV] default
    'quad_gradient': 16.0,       # [T/m] default
    'emitx_i': 1.e-6,            # [m*rad] Vorgabe emittance @ entrance
    'emity_i': 1.e-6,            # [m*rad] Vorgabe emittance @ entrance
    'betax_i': None,             # [m] Vorgabe twiss betax @ entrance
    'betay_i': None,             # [m] Vorgabe twiss betax @ entrance
    'alfax_i': None,             # Vorgabe twiss alphax @ entrance
    'alfay_i': None,             # Vorgabe twiss alphaxy @ entrance
    'dZ':      0.02,             # [m] Vorgabe longitudinal displacement dZ
    'dP/P':    0.02,             # [rad] Vorgabe relative impulse dP/P
    'dWf': False,                # acceleration on/off flag default
    'periodic': True,            # periodic lattice? default
    'verbose': 1,                # print flag (True) default
     }
CONF['wellenlänge']    = CONF['lichtgeschwindigkeit']/CONF['frequenz']
CONF['spalt_spannung'] = CONF['Ez_feld']*CONF['spalt_laenge']
SUMMARY = {}
class Particle(object):                           ## relativistic particle
    soll=None            #reference particle  (class member!)
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
        self.v          = self.beta*CONF['lichtgeschwindigkeit']    #   [m/s]
        self.brho       = 1.e6/CONF['lichtgeschwindigkeit']*self.p  #   [T*m]
        self.name       = name
    def incTK(self,deltaTK):
        self._set_self(self.tkin+deltaTK,self.e0,self.name)
    def out(self,tee=True):
        s=(u'          B*rho[Tm] Tk[MeV/c\u00B2] p[MeV/c]   gamma    beta     gamma*beta  E[MeV/c\u00B2]\n'+\
              '{:8s}{:8.4f}   {:8.4f}    {:8.4f} {:8.4f} {:8.4f}  {:8.4f}     {:8.4f}')\
            .format(self.name,self.brho,self.tkin,self.p,self.gamma,self.beta,self.gamma_beta,self.e)
        if tee:
            print(s)
        return s
    def TrTf(self,gap,fRF):  # transit-time-factor nach Panofsky (see Lapostolle CERN-97-09 pp.65)
        teta = 2.*pi*fRF*gap / (self.beta*CONF['lichtgeschwindigkeit'])
        teta = 0.5 * teta
        ttf = sin(teta)/teta
        return ttf

class Proton(Particle):                           ## proton
    def __init__(self,tkin=CONF['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass=CONF['proton_mass'],name='proton')

Particle.soll = Proton(CONF['injection_energy'])

def collect_summary():
	SUMMARY['frequency [Hz]'] = CONF['frequenz']
	SUMMARY['Quad gradient [T/m]'] = CONF['quad_gradient']
	SUMMARY['QF gradient [T/m]'] = CONF['quadf_gradient']
	SUMMARY['QD gradient [T/m]'] = CONF['quadd_gradient']
	SUMMARY['Quad pole length [m]'] = CONF['ql']
	SUMMARY['injection energy [MeV]'] = CONF['injection_energy']
	SUMMARY['emitx_i [rad*m]'] = CONF['emitx_i']
	SUMMARY['emity_i [rad*m]'] = CONF['emity_i']
	SUMMARY['sigx_i* [mm]'] = 1000.*sqrt(CONF['betax_i']*CONF['emitx_i'])  # enveloppe @ entrance
	SUMMARY['sigy_i* [mm]'] = 1000.*sqrt(CONF['betay_i']*CONF['emity_i'])
	SUMMARY['<impulse dP/P> [%]'] = CONF['dP/P'] * 1.e+2
	SUMMARY['sync. phase [deg]'] = CONF['soll_phase']
	SUMMARY['dZ [m]'] = CONF['dZ']
	SUMMARY['cavity gap length [m]'] = CONF['spalt_laenge']
	SUMMARY['cavity length [m]'] = CONF['cavity_laenge']
	SUMMARY['wavelength [m]'] = CONF['wellenlänge']
	SUMMARY['cavity gap voltage* [MV]'] = CONF['spalt_spannung']
	SUMMARY['acc. field Ez [MV/m]'] = CONF['Ez_feld']
	SUMMARY['lattice version'] = CONF['lattice_version']
	#...........*...........*...........*...........*...........*...........*
	SUMMARY['QF pole strength* [T]'] = CONF['quadf_gradient'] * CONF['ql']
	SUMMARY['QF current* [A/winding]'] = (CONF['quadf_gradient'] * (CONF['ql']*1000.)**2 )/2.52/CONF['n_coil']
	SUMMARY['QF power estimate* [W]'] = 0.0115 *SUMMARY['QF current* [A/winding]']**2  # R=0.0115 Ohms
	SUMMARY['QF coil [windings]'] = CONF['n_coil']
	SUMMARY['<energy dW/W> max* [%]'] = wakzp = Wakzeptanz(    # energy acceptance in %
		CONF['Ez_feld'],
		Particle.soll.TrTf(CONF['spalt_laenge'],CONF['frequenz']),
		CONF['soll_phase'],
		CONF['wellenlänge'],
		Particle.soll)*1.e+2
	SUMMARY['<impulse dP/P> max* [%]'] = 1./(1.+1./Particle.soll.gamma)*wakzp  # impule acceptanc in %
	SUMMARY['dphi* [deg]'] =degrees( 2*pi*CONF['frequenz']/CONF['lichtgeschwindigkeit']/Particle.soll.beta*CONF['dZ'])
	return

def Wakzeptanz(Ez,T,phis,lamb,particle):
    """
    Energieakzeptanz dW/W nach T.Wangler pp. 179
    Ez      [Mev/m] accelerating field gradient on gap axis
    T       transit time factor
    phis    [deg] sync. phase (-90 to 0)
    lamb    [m] rf wave length
    particle    sync. particle
    """
    gb = particle.gamma_beta
    m0 = particle.e0  # rest mass
    res = 2.*Ez*T*gb*gb*gb*lamb/pi/m0
    phsoll = radians(phis)
    res = res*(phsoll*cos(phsoll)-sin(phsoll))
    res = sqrt(res)/(particle.gamma-1)
    return res

class Electron(Particle):                         ## electron
    def __init__(self,tkin=CONF['injection_energy']):
        super(Electron,self).__init__(tkin=tkin,mass=CONF['electron_mass'],name='electron')

def k0(gradient=0.,tkin=0.):                  ## quad strength from B-field gradient & kin. energy
    """
    quad strength as function of kin. energy and gradient
    gradient: in [Tesla/m]
    tkin: kinetic energy in [MeV]
    (for protons only!)
    """
    if (tkin >= 0.):
        prot=Proton(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        kres = 1.e-6*CONF['lichtgeschwindigkeit']*gradient/(beta*gamma*e0)
        # print('k0= ',kres)
        return kres
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')

def scalek0(k0=0.,tki=0.,tkf=0.):             ## scale quad  strength with kin. energy
    """
    scale k0 for increase of kin. energy from
    tki to tkf  (for protons only!!)
    """
    pi  =Proton(tki)
    bi  =pi.beta
    gi  =pi.gamma
    pf  =Proton(tkf)
    bf  =pf.beta
    gf  =pf.gamma
    kf= k0 * (bi * gi) / (bf * gf)
    return kf

def dBdz(k0=0.,tkin=0.):                      ## B-field gradient from quad strength & kin. energy
    """
    calculate quad gradient for given quad strength k0
    and given kin. energy tkin
    (for protons only!!)
    """
    if (tkin >= 0.):
        prot=Proton(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        return k0*(beta*gamma*e0)/1.e-6*CONF['lichtgeschwindigkeit']
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')

def objprnt (what,text='========',filter=[]): ## helper to print objects as dictionary
    print('\n          ================= '+text+' =================')
    for k,v in sorted(what.__dict__.items()):
        if k in filter:
            continue
        print(k.rjust(30),':',v)
    return

def dictprnt(what,text='========',filter=[],njust=30): ## helper to print objects as dictionary
    print('\n          ================= '+text+' =================')
    for k,v in sorted(what.items()):
        if k in filter:
            continue
        print(k.rjust(njust),':',v)
    return

def printv(level,*args):                      ## multilevel printing using verbose flag
    verbose = CONF['verbose']
    if verbose >= level:
        print(*args)

def wille():
    return {
        'k_quad_f':1.2,
        'length_quad_f':0.2,
        'k_quad_d':1.2,
        'length_quad_d':0.4,
        'beding_radius':3.8197,
        'dipole_length':1.5,
        'drift_length':0.55
    }
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    Particle.soll = Particle()
    dictprnt(CONF,text='CONF')
    dictprnt(wille(),text='wille')
    print('\n'+Particle.soll.out(False))
    Proton(tkin=50.).out()
    Electron(tkin=50.).out()
    Particle.soll = Electron(50.)
    objprnt(Particle.soll,text='Particle.soll')
    Particle.soll = Proton(50.)
    objprnt(Particle.soll,text='Particle.soll')
