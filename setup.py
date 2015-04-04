# -*- coding: utf-8 -*-
import lattice as LAT
import elements as ELM

def values():
    return {
        'half_cell_length':2.6,
        'quad_length':1.243,
        'phase_advance_per_cell':108.,
        'bending_radius':279.38,
        'bending_magnet_length':2.486
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
def fodo():
    return {
        'k_quad_f':1.2,
        'length_quad_f':0.2,
        'k_quad_d':1.2,
        'length_quad_d':0.4,
        'beding_radius':3.8197,
        'dipole_length':1.5,
        'drift_length':0.55
        }

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
     md=ELM.D(ld)    
     ## lattice
     lattice=LAT.Lattice()
     lattice.add_element(mqf)
     lattice.add_element(md)
     lattice.add_element(mw)
     lattice.add_element(mb)
     lattice.add_element(mw)
     # lattice.add_element(mw1)
     # lattice.add_element(mb1)
     # lattice.add_element(mw1)
     lattice.add_element(md)
     lattice.add_element(mqd)
     lattice.add_element(md)
     lattice.add_element(mw)
     lattice.add_element(mb)
     lattice.add_element(mw)
     lattice.add_element(md)
     lattice.add_element(mqf)
     # lattice.out()
     top=LAT.Lattice()
     top.append(lattice)
     # top.append(lattice)
     # top.append(lattice)
     # top.append(lattice)
     # top.append(lattice)
     top.append(top)
     top.append(top)
     top.append(top)
     top.append(top)
     top.append(top)
     top.append(top)
     return top

if __name__ == '__main__':
     lat = make_lattice()
     lat.out()
     print(lat)
