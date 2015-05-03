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
def objprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
        print('\n========= '+text+' =================')
        for k,v in what.__dict__.items():
            if k in filter:
                continue
            print(k.rjust(30),':',v)
def dictprnt(what,text='========',filter={}):  ## helper to print objects as dictionary
        print('\n========= '+text+' =================')
        for k,v in what.items():
            if k in filter:
                continue
            print(k.rjust(30),':',v)
