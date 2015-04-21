# -*- coding: utf-8 -*-
from math import pi,sqrt
Phys = {
    'lichtgeschwindigkeit': 299792458.,    #  m/s
    'elementarladung': 1.602176565e-19,  # coulomb
    'proton_mass': 938.272,  # MeV/c**2
    'spalt_spannung': 3.0,   # MV
    'spalt_laenge':0.04,     # m
    'transit_time': 0.5,
    'soll_phase': -50.0,     # deg
    'frequenz': 800.0,       # MHz
    'kinetic_energy': 50.,   #MeV
    'quad_gradient': 1.0,    # T/m
    'radians': pi/180.,      # rad/deg
    'degrees': 180./pi,      # deg/rad
     }
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
class Beam():
    def __init__(self,tkin=0.):
        self.e0   = Phys['proton_mass']   # proton rest mass [MeV/c**2]
        self.tkin = tkin                     # proton kinetic energy [MeV]
        self.e    = self.e0+self.tkin        # proton total energy [MeV]
        self.gamma= self.e/self.e0
        self.beta = sqrt(1.-1./(self.gamma*self.gamma))
        self.v    = self.beta*Phys['lichtgeschwindigkeit']
        self.name = 'proton'
    def out(self):
        print('{:s}:  T-kin[MeV]={:.3f} gamma {:.3f} beta {:.3f} velocity[m/s] {:.6g} E[MeV] {:.3f} '
            .format(self.name,self.tkin,self.gamma,self.beta,self.v,self.e))
def test1():
    print('\ntest: Beam class')
    for key,value in Phys.items():
        print('{}=  {:.4g}'.format(key.rjust(20),value))
    # Beam class
    print()
    Beam(0.).out()
    Beam(50.).out()
    Beam(200.).out()
    Beam(1.e6).out()
    Beam(1.e9).out()
if __name__ == '__main__':
    test1()