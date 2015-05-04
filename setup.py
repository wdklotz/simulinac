 # -*- coding: utf-8 -*-
from math import pi,sqrt

Phys = {
    'lichtgeschwindigkeit': 299792458.,    # [m/s]
    'elementarladung': 1.602176565e-19,    # [coulomb]
    'proton_mass': 938.272,  # [MeV/c**2]
    'spalt_spannung': 3.5,   # [MV]
    'spalt_laenge': 0.08,    # [m]
    'soll_phase': -45.0,     # [deg]
    'frequenz': 800.0,       # [MHz]
    'injection_energy': 50., # [MeV]
    'quad_gradient': 1.0,    # [T/m]
    'radians': pi/180.,      # [rad/deg]
    'degrees': 180./pi,      # [deg/rad]
    'emitx_i':5.e-6,         # [m*rad] emttance @ entrance
    'emity_i':5.e-6,         # [m*rad] emttance @ entrance
    'sigx_i': 5.e-3,         # [m] one sigma transverse beam size
    'sigy_i': 2.5e-3,        # [m] one sigma transverse beam size
    'dZ': 1.8e-2,            # [m] longitudinal displacement dZ
    'dP/P': 5.e-2,           # [rad] relative impulse dP/P
     }
Phys['wellenlänge']=1.e-6*Phys['lichtgeschwindigkeit']/Phys['frequenz']
class Beam():   ## relativistic protons
    soll=None   ## the synchronous reference particle  (class member!)
    def __init__(self,tkin=0.):
        self._set_self(tkin)
    def _set_self(self,tkin):
        self.tkin = tkin                     # proton kinetic energy [MeV]
        self.e0   = Phys['proton_mass']      # proton rest mass [MeV/c**2]
        self.e    = self.e0+self.tkin        # proton total energy [MeV]
        self.gamma= self.e/self.e0
        self.beta = sqrt(1.-1./(self.gamma*self.gamma))
        self.v    = self.beta*Phys['lichtgeschwindigkeit']
        self.name = 'proton'
    def incTK(self,deltaTK):
        self._set_self(self.tkin+deltaTK)
    def out(self):
        print('{:s}:  T-kin[MeV]={:.3f} gamma {:.3f} beta {:.3f} velocity[m/s] {:.6g} E[MeV] {:.3f} '
            .format(self.name,self.tkin,self.gamma,self.beta,self.v,self.e))        
Beam.soll = Beam(Phys['injection_energy']) # the synchronous reference particle  (class member!)
def k0(gradient=0.,tkin=0.):       ## quad strength from B-field gradient & kin. energy
    """
    quad strength as function of kin. energy and gradient
    gradient: in [Tesla/m]
    tkin: kinetic energy in [MeV]
    """
    if (tkin >= 0.):
        prot=Beam(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*Phys['lichtgeschwindigkeit']
        kres = factor*gradient/(beta*gamma*e0) 
        # print('k0= ',kres)
        return kres
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
def scalek0(k0=0.,tki=0.,tkf=0.):  ## scale quad  strength with kin. energy
    """
    scale k0 for increase of kin. energy from
    tki to tkf
    """
    pi  =Beam(tki)
    bi  =pi.beta
    gi  =pi.gamma
    pf  =Beam(tkf)
    bf  =pf.beta
    gf  =pf.gamma
    kf= k0 * (bi * gi) / (bf * gf)
    return kf
def dBdz(k0=0.,tkin=0.):           ## B-field gradient from quad strength & kin. energy
    """
    calculate quad gradient for given quad strength k0
    and given kin. energy tkin
    """
    if (tkin >= 0.):
        prot=Beam(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*Phys['lichtgeschwindigkeit']
        return k0*(beta*gamma*e0)/factor           
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
def objprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
    print('\n========= '+text+' =================')
    for k,v in what.__dict__.items():
        if k in filter:
            continue
        print(k.rjust(30),':',v)
    return
def dictprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
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
    dictprnt(Phys,'Phys')
    objprnt(Beam.soll,'Beam.soll')
    dictprnt(wille(),'wille')