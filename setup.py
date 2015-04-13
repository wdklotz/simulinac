# -*- coding: utf-8 -*-
import lattice as LAT
import elements as ELM
from math import pi, sqrt

physics = {
    'lichtgeschwindigkeit': 299792458.,    #  m/s
    # 'elementarladung': 1.602176565e-19,  # coulomb
    'proton_mass': 938.272,  # MeV/c**2
    'spalt_spannung': 5.0,   # MV
    'transit_time': 0.5,
    'soll_phase': -60.0,     # deg
    'frequenz': 800.0,       # MHz
    'kinetic_energy': 50.,   #MeV
    'quad_gradient': 1.0,    # T/m
    'quad_gradient': 1.0,    # T/m
    'radians': pi/180.,      # rad/deg
    'degrees': 180./pi,      # deg/rad
     }
def ex_wille():
    return {
        'k_quad_f':1.2,
        'length_quad_f':0.2,
        'k_quad_d':1.2,
        'length_quad_d':0.4,
        'beding_radius':3.8197,
        'dipole_length':1.5,
        'drift_length':0.55
    }
class Proton():
    def __init__(self,tkin=0.):
        self.e0   = physics['proton_mass']   #  proton rest mass [MeV/c**2]
        self.tkin = tkin      # proton kinetic energy
        self.e    = self.e0+self.tkin  # proton total energy
        self.gamma= self.e/self.e0
        self.beta = sqrt(1.-1./(self.gamma))
        self.v    = self.beta*physics['lichtgeschwindigkeit']
        self.name = 'proton'
    def out(self):
        print('{:s}:  T-kin[MeV]={:.3f} gamma {:.3f} beta {:.3f} velocity[m/s] {:.6g} E[MeV] {:.3f} '
            .format(self.name,self.tkin,self.gamma,self.beta,self.v,self.e))
def k0p(gradient=0.,tkin=0.):
    """
    quad strength as function of kin. energy and gradient
    particle: proton
    gradient: in [Tesla/m]
    tkin: kinetic energy in [MeV]
    """
    if (tkin >= 0.):
        prot=Proton(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*physics['lichtgeschwindigkeit']
        return factor*gradient/(beta*gamma*e0)             
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
def bgrad(k0=0.,tkin=0.):
    if (tkin >= 0.):
        prot=Proton(tkin)
        beta=prot.beta
        e0=prot.e0
        gamma=prot.gamma
        factor=1.e-6*physics['lichtgeschwindigkeit']
        return k0*(beta*gamma*e0)/factor           
    else:
        raise RuntimeError('setup.k0(): negative kinetic energy?')
    
def make_lattice():
     print("K.Wille's Beispiel auf pp. 112-113")
     kqf=  ex_wille()['k_quad_f']
     lqf=  ex_wille()['length_quad_f']
     kqd=  ex_wille()['k_quad_d']
     lqd=  ex_wille()['length_quad_d']
     rhob= ex_wille()['beding_radius'] 
     lb=   ex_wille()['dipole_length']
     ld=   ex_wille()['drift_length']
     ## elements
     mqf=ELM.QF(kqf,lqf,'QF')
     mqd=ELM.QD(kqd,lqd,'QD')
     mb=ELM.SD(rhob,lb,'B')
     mb1=ELM.SD(rhob,lb*0.5,'B1')  ## 1/2 sector dip.
     mw=ELM.WD(lb,rhob)
     mw1=ELM.WD(lb*0.5,rhob)
     mbr=ELM.RD(rhob,lb)
     md=ELM.D(ld)    
     ## lattice
     lattice=LAT.Lattice()
     lattice.add_element(mqf)
     lattice.add_element(md)
     # lattice.add_element(mw)
     # lattice.add_element(mb)
     # lattice.add_element(mw)
     lattice.add_element(mbr)
     # lattice.add_element(mw1)
     # lattice.add_element(mb1)
     # lattice.add_element(mw1)
     lattice.add_element(md)
     lattice.add_element(mqd)
     lattice.add_element(md)
     # lattice.add_element(mw)
     # lattice.add_element(mb)
     # lattice.add_element(mw)
     lattice.add_element(mbr)
     lattice.add_element(md)
     lattice.add_element(mqf)
     # lattice.out()
     top=LAT.Lattice()
     top.append(lattice)
     top.append(lattice)
     top.append(lattice)
     top.append(lattice)
     top.append(lattice)
     top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     return top
def test0():
    lat = make_lattice()
    lat.out()
    print(lat)
def test1():
    print('\ntest: Proton class')
    for key,value in physics.items():
        print('{}=  {:.4f}'.format(key.rjust(20),value))
    # Proton class
    print()
    Proton(0.).out()
    Proton(50.).out()
    Proton(200.).out()
    Proton(1.e6).out()
    Proton(1.e9).out()
def test2():
    print('\ntest: quad k-faktor')
    Bgrad=physics['quad_gradient']   # [T/m] gradient
    tk=physics['kinetic_energy']      # [MeV]  kin. energy
    kq=k0p(gradient=Bgrad,tkin=tk)   # quad strength [1/m**2]
    len=0.2   # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]
    print('\nproton {:.3f}[MeV] ~ beta {:.3f} in quadrupole:'.format(tk,Proton(tk).beta))
    print('k [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(Bgrad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))
if __name__ == '__main__':
    test0()
    test1()
    test2()
