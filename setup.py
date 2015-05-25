 # -*- coding: utf-8 -*-
from math import pi,sqrt,sin

CONF = {                                      ## CONFIG constants and setup ...
    'lichtgeschwindigkeit': 299792458.,    # [m/s]
    'elementarladung': 1.602176565e-19,    # [coulomb]
    'proton_mass': 938.272,      # [MeV/c**2]
    'electron_mass': 0.5109989,  # [MeV/c**2]
    'Ez_feld': 4.,               # [MV/m]
    'spalt_laenge': 0.04,        # [m]
    'cavity_laenge': 0.08,       # [m]
    'soll_phase': -45.0,         # [deg]
    'frequenz': 800.e6,          # [Hz]
    'injection_energy': 50.,     # [MeV]
    'quad_gradient': 7.5,        # [T/m]
    'emitx_i':5.e-6,             # [m*rad] emttance @ entrance
    'emity_i':5.e-6,             # [m*rad] emttance @ entrance
    'sigx_i': 5.e-3,             # [m] one sigma transverse beam size
    'sigy_i': 2.5e-3,            # [m] one sigma transverse beam size
    'dZ': 1.8e-2,                # [m] longitudinal displacement dZ
    'dP/P': 2.e-2,               # [rad] relative impulse dP/P
    'dWf': 1.0,                  # acceleration on/off flag
    'periodic': True,            # periodic lattice?
     }
CONF['wellenlänge']   =CONF['lichtgeschwindigkeit']/CONF['frequenz']
CONF['spalt_spannung']=CONF['Ez_feld']*CONF['spalt_laenge']
class Beam(object):                           ## relativistic particle beam
    soll=None   ## the synchronous reference particle  (class member!)
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
class Proton(Beam):                           ## proton
    def __init__(self,tkin=CONF['injection_energy']):
        super(Proton,self).__init__(tkin=tkin,mass=CONF['proton_mass'],name='proton')
Beam.soll             = Proton(CONF['injection_energy'])
class Electron(Beam):                         ## electron
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
    print('\n========= '+text+' =================')
    for k,v in what.__dict__.items():
        if k in filter:
            continue
        print(k.rjust(30),':',v)
    return
def dictprnt(what,text='========',filter=[]): ## helper to print objects as dictionary
    print('\n========= '+text+' =================')
    for k,v in what.items():
        if k in filter:
            continue
        print(k.rjust(30),':',v)
    return
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
    Beam.soll = Beam()
    dictprnt(CONF,text='CONF')
    dictprnt(wille(),text='wille')
    print('\n'+Beam.soll.out(False))
    Proton(tkin=50.).out()
    Electron(tkin=50.).out()
    Beam.soll = Electron(50.)
    objprnt(Beam.soll,text='Beam.soll')
    Beam.soll = Proton(50.)
    objprnt(Beam.soll,text='Beam.soll')
