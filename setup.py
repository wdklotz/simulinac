# -*- coding: utf-8 -*-
import lattice as LAT
import elements as ELM
from math import pi, sqrt

physics = {
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
        self.e0   = physics['proton_mass']   # proton rest mass [MeV/c**2]
        self.tkin = tkin                     # proton kinetic energy [MeV]
        self.e    = self.e0+self.tkin        # proton total energy [MeV]
        self.gamma= self.e/self.e0
        self.beta = sqrt(1.-1./(self.gamma*self.gamma))
        self.v    = self.beta*physics['lichtgeschwindigkeit']
        self.name = 'proton'
    def out(self):
        print('{:s}:  T-kin[MeV]={:.3f} gamma {:.3f} beta {:.3f} velocity[m/s] {:.6g} E[MeV] {:.3f} '
            .format(self.name,self.tkin,self.gamma,self.beta,self.v,self.e))
def k0(gradient=0.,tkin=0.):
    """
    quad strength as function of kin. energy and gradient
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
def k0scaled(k0=0.,tkin0=0.,tkin1=0.):
    """
    scale k0 for increase of kin. energy from
    tkin0 to tkin1
    """
    prot0  =Proton(tkin0)
    prot1  =Proton(tkin1)
    beta0  =prot0.beta
    gamma0 =prot0.gamma
    beta1  =prot1.beta
    gamma1 =prot1.gamma
    k1= k0 * (beta0 * gamma0) / (beta1 * gamma1)
    return k1
def QDscaled(quad0,tkin0=0.,tkin1=0.):
    k0   =quad0.k0
    len  =quad0.length
    label=quad0.label
    k1=k0scaled(k0,tkin0,tkin1)
    if isinstance(quad0,ELM.QF) and (isinstance(quad0,ELM.QD)==False):
        quad1=ELM.QF(k0=k1,length=len,label=label)
    elif isinstance(quad0,ELM.QD):
        quad1=ELM.QD(k0=k1,length=len,label=label)
    else:
        raise RuntimeError('QF._mx: neither QF nor QD! should never happen!')
    return quad1
def CAVscaled(cavity,tkin=0.):
    u0=cavity.u0
    phiSoll=cavity.phis
    fRF=cavity.freq
    label=cavity.label
    cav = ELM.CAV(U0=u0,PhiSoll=phiSoll,Tkin=tkin,fRF=fRF,label=label)
    return cav
def dBdz_p(k0=0.,tkin=0.):
    """
    calculate quad gradient for given quad strength k0
    and given kin. energy tkin
    """
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
    kq=k0(gradient=Bgrad,tkin=tk)   # quad strength [1/m**2]
    len=0.4   # quad len [m]
    focal = kq*len
    focal=1./focal  # focal len [m]

    print('\nproton {:.3f}[MeV] ~ beta {:.3f} in quadrupole:'.format(tk,Proton(tk).beta))
    print('k [1/m**2]\t{:3f}'.format(kq))
    print('dB/dz[T/m]\t{:.3f}'.format(Bgrad))
    print('len[m]\t\t{:.3f}'.format(len))
    print('focal len[m]\t{:.3f}'.format(focal))
    
    Bgrad=dBdz_p(kq,tk)
    print('dB/dz[T/m]\t{:.3f}'.format(Bgrad))
    mqf=ELM.QF(kq,len)
    mqd=ELM.QD(kq,len)
    cavity=ELM.CAV(
        # U0=physics['spalt_spannung'],
        U0=3.0,
        # TrTF=physics['transit_time'],
        PhiSoll=physics['soll_phase']*physics['radians'],
        Tkin=physics['kinetic_energy'],
        fRF=physics['frequenz'],
        label='gap')
    for dt in [10.,50.,150.]:
        tks=tk+dt
        k_scaled = k0scaled(kq,tk,tks)
        print('k[{} MeV] {:.3f} k[{} MeV] {:.3f}'.format(tk,kq,tks,k_scaled))
        QDscaled(mqf,tk,tks).out()
        QDscaled(mqd,tk,tks).out()
        CAVscaled(cavity,tks).out()
if __name__ == '__main__':
    test0()
    test1()
    test2()
