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
############################################################################